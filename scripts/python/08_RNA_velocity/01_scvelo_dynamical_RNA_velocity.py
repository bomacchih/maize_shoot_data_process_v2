#!/usr/bin/env python3
"""Dynamical RNA-velocity analysis for five maize shoot domains.

Expected project structure
--------------------------
data/
├── metadata/metadata.csv
└── processed/
    ├── XGE20_UL01_B73V5S210_short_MTCL.loom
    ├── XGE20_UL02_B73V5S210_short_MTCL.loom
    ├── XGE20_UL04_B73V5S210_short_MTCL.loom
    ├── XGE21_VR01_B73V5S210_short_MTCL.loom
    └── ... one loom file for each of the 14 samples

The loom files must have been generated from the corresponding Space Ranger
outputs with Velocyto and the Zm-B73-REFERENCE-NAM-5.0 GTF. This program does
not run shell commands or Velocyto. It performs the Python/scVelo portion only.

Only spots labeled SAM, P1_P2, P3, P4, or P5 in metadata.csv are retained.
Overlapping-domain spots are therefore excluded by the curated metadata table.
"""

from __future__ import annotations

import argparse
import importlib.metadata
import json
import logging
import re
import sys
from pathlib import Path

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import scvelo as scv
from scipy import sparse


SAMPLE_TO_SUFFIX = {
    "UL01": 1,
    "UL02": 2,
    "UL04": 3,
    "VR01": 4,
    "VR02": 5,
    "VR03": 6,
    "VR04": 7,
    "DQ01": 8,
    "DQ02": 9,
    "DQ03": 10,
    "DQ04": 11,
    "DQ06": 12,
    "DQ07": 13,
    "DQ08": 14,
}

DOMAIN_ORDER = ["SAM", "P1_P2", "P3", "P4", "P5"]
DOMAIN_COLORS = {
    "SAM": "#E41A1C",
    "P1_P2": "#F8766D",
    "P3": "#2CA25F",
    "P4": "#19BFC4",
    "P5": "#2C6DB2",
}


def project_root_from_script() -> Path:
    """Return repository root for scripts/python/08_RNA_velocity/script.py."""
    return Path(__file__).resolve().parents[3]


def parse_arguments() -> argparse.Namespace:
    root = project_root_from_script()
    parser = argparse.ArgumentParser(
        description="Run dynamical scVelo analysis for maize developing leaves."
    )
    parser.add_argument("--project-root", type=Path, default=root)
    parser.add_argument(
        "--raw-dir",
        type=Path,
        default=None,
        help=(
            "Directory containing the loom files. Defaults to data/raw/loom, "
            "where the Zenodo record should be downloaded. "
            "Legacy per-sample subdirectories are also supported."
        ),
    )
    parser.add_argument("--metadata", type=Path, default=None)
    parser.add_argument("--processed-dir", type=Path, default=None)
    parser.add_argument("--figure-dir", type=Path, default=None)
    parser.add_argument("--table-dir", type=Path, default=None)
    parser.add_argument(
        "--min-shared-counts",
        type=int,
        default=20,
        help="Minimum combined spliced/unspliced counts used by scVelo gene filtering.",
    )
    parser.add_argument(
        "--n-top-genes",
        type=int,
        default=2000,
        help="Number of highly variable genes retained for velocity modeling.",
    )
    parser.add_argument(
        "--n-pcs",
        type=int,
        default=30,
        help="Number of PCs used for the neighborhood graph and moments.",
    )
    parser.add_argument(
        "--n-neighbors",
        type=int,
        default=30,
        help="Number of neighboring spots used to calculate moments.",
    )
    parser.add_argument(
        "--max-iter",
        type=int,
        default=20,
        help="Maximum EM iterations used by recover_dynamics().",
    )
    parser.add_argument(
        "--n-jobs",
        type=int,
        default=1,
        help="Parallel jobs for recover_dynamics(); use 1 for maximum portability.",
    )
    parser.add_argument(
        "--use-harmony",
        action="store_true",
        help=(
            "Optionally correct the PCA representation by sample_id before moments. "
            "Off by default because Harmony was not specified for RNA velocity in "
            "the manuscript method."
        ),
    )
    parser.add_argument(
        "--diagnostic-unspliced-threshold",
        type=float,
        default=0.01,
        help=(
            "Report sample/domain groups with median unspliced fraction below this "
            "value. This diagnostic does not remove spots."
        ),
    )
    parser.add_argument(
        "--validate-files-only",
        action="store_true",
        help=(
            "Validate the 14 loom filenames and metadata coverage, then exit "
            "without loading loom matrices or running scVelo."
        ),
    )
    return parser.parse_args()


def configure_paths(args: argparse.Namespace) -> argparse.Namespace:
    root = args.project_root.resolve()
    args.project_root = root
    args.raw_dir = (args.raw_dir or root / "data" / "raw" / "loom").resolve()
    args.metadata = (
        args.metadata or root / "data" / "metadata" / "metadata.csv"
    ).resolve()
    args.processed_dir = (
        args.processed_dir or root / "data" / "processed" / "RNA_velocity"
    ).resolve()
    args.figure_dir = (
        args.figure_dir or root / "results" / "figures" / "08_RNA_velocity"
    ).resolve()
    args.table_dir = (
        args.table_dir or root / "results" / "tables" / "08_RNA_velocity"
    ).resolve()

    for directory in (args.processed_dir, args.figure_dir, args.table_dir):
        directory.mkdir(parents=True, exist_ok=True)
    return args


def configure_logging(table_dir: Path) -> None:
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
        handlers=[
            logging.StreamHandler(sys.stdout),
            logging.FileHandler(table_dir / "RNA_velocity_run.log", mode="w"),
        ],
    )


def find_one_loom(raw_dir: Path, sample_id: str) -> Path:
    """Find exactly one sample loom in the flat or legacy directory layout."""
    if not raw_dir.is_dir():
        raise FileNotFoundError(f"Loom directory not found: {raw_dir}")

    candidates: list[Path] = []

    # Canonical layout: all loom files are directly under data/raw/loom/.
    candidates.extend(
        path for path in raw_dir.glob("*.loom") if sample_id.lower() in path.name.lower()
    )

    # Backward compatibility: allow data/raw/<sample>/.../*.loom or another
    # nested archive layout when --raw-dir points to a parent directory.
    if not candidates:
        candidates.extend(
            path
            for path in raw_dir.rglob("*.loom")
            if sample_id.lower() in path.name.lower()
        )

    candidates = sorted({path.resolve() for path in candidates})
    if len(candidates) == 0:
        raise FileNotFoundError(
            f"No loom file found for {sample_id} under {raw_dir}. "
            "Expected a filename containing the sample ID, such as "
            f"data/raw/loom/XGE20_{sample_id}_B73V5S210_short_MTCL.loom."
        )
    if len(candidates) > 1:
        formatted = "\n  ".join(str(path) for path in candidates)
        raise RuntimeError(
            f"More than one loom file was found for {sample_id}:\n  {formatted}\n"
            "Keep one loom per sample or place each sample in its own directory."
        )
    return candidates[0]


def canonical_barcode(raw_barcode: str, sample_number: int) -> str:
    """Convert a Velocyto barcode to the Seurat-v5 metadata barcode format."""
    barcode = str(raw_barcode).split(":")[-1].strip()

    # Velocyto loom files commonly end barcodes with "x" rather than "-1".
    barcode = re.sub(r"x$", "", barcode)

    # Remove a previously appended merged-object suffix if present.
    barcode = re.sub(r"_1_[0-9]+$", "", barcode)

    if not barcode.endswith("-1"):
        barcode = f"{barcode}-1"
    return f"{barcode}_1_{sample_number}"


def set_stable_gene_ids(loom: ad.AnnData, sample_id: str) -> None:
    """Prefer stable gene accessions when Velocyto stored symbols as var_names."""
    for column in ("Accession", "gene_id", "Gene"):
        if column not in loom.var.columns:
            continue
        values = loom.var[column].astype(str)
        if values.notna().all() and values.is_unique:
            loom.var_names = values
            break

    if loom.var_names.has_duplicates:
        logging.warning(
            "%s contains duplicated gene identifiers; making them unique.", sample_id
        )
        loom.var_names_make_unique()


def load_metadata(metadata_path: Path) -> pd.DataFrame:
    if not metadata_path.is_file():
        raise FileNotFoundError(f"Metadata file not found: {metadata_path}")

    metadata = pd.read_csv(metadata_path, dtype=str)
    required = {"Barcode", "sample_id", "domains"}
    missing = required.difference(metadata.columns)
    if missing:
        raise ValueError(
            "Metadata is missing required column(s): " + ", ".join(sorted(missing))
        )
    if metadata["Barcode"].duplicated().any():
        duplicated = metadata.loc[metadata["Barcode"].duplicated(), "Barcode"].head()
        raise ValueError(
            "metadata.csv contains duplicated Barcode values, for example: "
            + ", ".join(duplicated)
        )

    metadata = metadata.loc[metadata["domains"].isin(DOMAIN_ORDER)].copy()
    metadata = metadata.loc[metadata["sample_id"].isin(SAMPLE_TO_SUFFIX)].copy()
    metadata["domains"] = pd.Categorical(
        metadata["domains"], categories=DOMAIN_ORDER, ordered=True
    )
    metadata = metadata.set_index("Barcode", drop=False)

    if metadata.empty:
        raise ValueError("No SAM/P1_P2/P3/P4/P5 spots remained after metadata filtering.")
    return metadata


def validate_input_files(raw_dir: Path, metadata: pd.DataFrame) -> None:
    """Validate the inexpensive parts of the input contract before analysis."""
    for sample_id in SAMPLE_TO_SUFFIX:
        loom_path = find_one_loom(raw_dir, sample_id)
        if loom_path.stat().st_size == 0:
            raise ValueError(f"Loom file is empty: {loom_path}")

        target_spots = int((metadata["sample_id"] == sample_id).sum())
        if target_spots == 0:
            raise ValueError(
                f"metadata.csv has no retained SAM/P1_P2/P3/P4/P5 spots for {sample_id}."
            )
        logging.info(
            "Validated %s: %s (%d bytes; %d target metadata spots)",
            sample_id,
            loom_path.name,
            loom_path.stat().st_size,
            target_spots,
        )


def load_and_subset_looms(
    raw_dir: Path, metadata: pd.DataFrame
) -> tuple[ad.AnnData, pd.DataFrame]:
    looms: list[ad.AnnData] = []
    matching_rows: list[dict[str, object]] = []

    for sample_id, sample_number in SAMPLE_TO_SUFFIX.items():
        loom_path = find_one_loom(raw_dir, sample_id)
        logging.info("Loading %s: %s", sample_id, loom_path)
        # scVelo 0.3.x no longer exposes scv.read(). Scanpy's current loom
        # reader returns an AnnData object and preserves sparse count layers.
        loom = sc.read_loom(str(loom_path), sparse=True)

        required_layers = {"spliced", "unspliced"}
        missing_layers = required_layers.difference(loom.layers.keys())
        if missing_layers:
            raise ValueError(
                f"{loom_path} is missing layer(s): {', '.join(sorted(missing_layers))}"
            )

        set_stable_gene_ids(loom, sample_id)
        original_spots = loom.n_obs
        loom.obs_names = [
            canonical_barcode(barcode, sample_number) for barcode in loom.obs_names
        ]
        if loom.obs_names.has_duplicates:
            raise ValueError(f"Canonical barcodes are duplicated in {loom_path}")

        expected = metadata.index[metadata["sample_id"] == sample_id]
        matched = expected.intersection(loom.obs_names, sort=False)
        missing = expected.difference(loom.obs_names, sort=False)

        matching_rows.append(
            {
                "sample_id": sample_id,
                "loom_file": str(loom_path),
                "spots_in_loom": original_spots,
                "target_metadata_spots": len(expected),
                "matched_spots": len(matched),
                "missing_target_spots": len(missing),
            }
        )

        if len(matched) == 0:
            raise ValueError(
                f"No {sample_id} loom barcodes matched the curated metadata. "
                "Check the sample-to-suffix mapping and loom barcode format."
            )
        if len(missing) > 0:
            logging.warning(
                "%s: %d of %d target metadata spots were absent from the loom file.",
                sample_id,
                len(missing),
                len(expected),
            )

        loom = loom[matched, :].copy()
        loom.obs["loom_sample_id"] = sample_id
        looms.append(loom)

    combined = ad.concat(
        looms,
        axis=0,
        join="inner",
        merge="same",
        uns_merge="same",
        index_unique=None,
    )
    if combined.obs_names.has_duplicates:
        raise ValueError("Duplicated spot barcodes were found after loom concatenation.")

    ordered_spots = metadata.index.intersection(combined.obs_names, sort=False)
    combined = combined[ordered_spots, :].copy()
    combined.obs = metadata.loc[ordered_spots].copy()
    combined.obs["domains"] = pd.Categorical(
        combined.obs["domains"], categories=DOMAIN_ORDER, ordered=True
    )

    return combined, pd.DataFrame(matching_rows)


def row_sum(matrix: object) -> np.ndarray:
    values = matrix.sum(axis=1)
    if sparse.issparse(values):
        values = values.A
    return np.asarray(values).ravel()


def add_unspliced_qc(adata: ad.AnnData) -> None:
    spliced = row_sum(adata.layers["spliced"])
    unspliced = row_sum(adata.layers["unspliced"])
    total = spliced + unspliced
    fraction = np.divide(
        unspliced,
        total,
        out=np.zeros_like(unspliced, dtype=float),
        where=total > 0,
    )
    adata.obs["spliced_counts"] = spliced
    adata.obs["unspliced_counts"] = unspliced
    adata.obs["unspliced_fraction"] = fraction


def export_qc(
    adata: ad.AnnData,
    matching: pd.DataFrame,
    table_dir: Path,
    diagnostic_threshold: float,
) -> None:
    matching.to_csv(table_dir / "loom_metadata_barcode_matching.csv", index=False)

    summary = (
        adata.obs.groupby(["sample_id", "domains"], observed=True)
        .agg(
            n_spots=("Barcode", "size"),
            median_spliced_counts=("spliced_counts", "median"),
            median_unspliced_counts=("unspliced_counts", "median"),
            median_unspliced_fraction=("unspliced_fraction", "median"),
        )
        .reset_index()
    )
    summary["low_unspliced_diagnostic"] = (
        summary["median_unspliced_fraction"] < diagnostic_threshold
    )
    summary.to_csv(table_dir / "spliced_unspliced_QC_by_sample_domain.csv", index=False)

    low_groups = summary.loc[summary["low_unspliced_diagnostic"]]
    if not low_groups.empty:
        logging.warning(
            "%d sample-domain group(s) have median unspliced fraction below %.3f. "
            "They were reported but not automatically removed.",
            len(low_groups),
            diagnostic_threshold,
        )


def use_existing_umap_if_available(adata: ad.AnnData) -> bool:
    candidates = [
        ("UMAP_1", "UMAP_2"),
        ("umap_1", "umap_2"),
        ("harmonyUMAP_1", "harmonyUMAP_2"),
    ]
    for x_column, y_column in candidates:
        if x_column not in adata.obs or y_column not in adata.obs:
            continue
        coordinates = adata.obs[[x_column, y_column]].apply(pd.to_numeric, errors="coerce")
        if coordinates.notna().all().all():
            adata.obsm["X_umap"] = coordinates.to_numpy(dtype=float)
            logging.info("Using metadata UMAP coordinates: %s, %s", x_column, y_column)
            return True
    return False


def preprocess_and_run_velocity(adata: ad.AnnData, args: argparse.Namespace) -> None:
    scv.pp.filter_and_normalize(
        adata,
        min_shared_counts=args.min_shared_counts,
        n_top_genes=args.n_top_genes,
    )

    n_pcs = min(args.n_pcs, adata.n_obs - 1, adata.n_vars - 1)
    if n_pcs < 2:
        raise ValueError("Too few spots or genes remain for PCA and velocity analysis.")

    # Compute PCA explicitly so the representation used for moments is recorded.
    sc.tl.pca(
        adata,
        n_comps=n_pcs,
        use_highly_variable=True,
        svd_solver="arpack",
        random_state=0,
    )
    representation = "X_pca"

    if args.use_harmony:
        try:
            import scanpy.external as sce

            sce.pp.harmony_integrate(
                adata,
                key="sample_id",
                basis="X_pca",
                adjusted_basis="X_pca_harmony",
                random_state=0,
            )
        except ImportError as error:
            raise ImportError(
                "--use-harmony requires the optional harmonypy package."
            ) from error
        representation = "X_pca_harmony"

    scv.pp.moments(
        adata,
        n_neighbors=min(args.n_neighbors, adata.n_obs - 1),
        n_pcs=n_pcs if representation == "X_pca" else None,
        use_rep=representation,
    )

    if "X_umap" not in adata.obsm:
        sc.tl.umap(adata, random_state=0)

    # Dynamical mode requires kinetic-parameter recovery before velocity.
    scv.tl.recover_dynamics(
        adata,
        n_top_genes=min(args.n_top_genes, adata.n_vars),
        max_iter=args.max_iter,
        n_jobs=args.n_jobs,
    )
    scv.tl.velocity(adata, mode="dynamical")
    scv.tl.velocity_graph(adata, n_jobs=args.n_jobs)
    scv.tl.velocity_confidence(adata)


def save_current_figure(path: Path) -> None:
    figure = plt.gcf()
    figure.savefig(path, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(figure)


def generate_figures(adata: ad.AnnData, figure_dir: Path) -> None:
    palette = [DOMAIN_COLORS[domain] for domain in DOMAIN_ORDER]

    sc.pl.umap(
        adata,
        color="domains",
        palette=palette,
        frameon=False,
        show=False,
        title="Developing shoot domains",
    )
    save_current_figure(figure_dir / "Figure_10A_domains_UMAP.png")

    scv.pl.proportions(adata, groupby="domains", show=False)
    save_current_figure(figure_dir / "RNA_velocity_spliced_unspliced_proportions.png")

    scv.pl.velocity_embedding_grid(
        adata,
        basis="umap",
        color="domains",
        palette=palette,
        density=1.0,
        arrow_length=3,
        arrow_size=2,
        frameon=False,
        title="",
        show=False,
    )
    save_current_figure(figure_dir / "Figure_10B_scVelo_dynamical_velocity_grid.png")

    scv.pl.velocity_embedding_stream(
        adata,
        basis="umap",
        color="domains",
        palette=palette,
        density=2,
        max_length=1,
        frameon=False,
        title="",
        show=False,
    )
    save_current_figure(figure_dir / "scVelo_dynamical_velocity_stream.png")


def package_versions() -> dict[str, str]:
    packages = [
        "python",
        "anndata",
        "scanpy",
        "scvelo",
        "numpy",
        "pandas",
        "scipy",
        "matplotlib",
        "loompy",
    ]
    versions: dict[str, str] = {"python": sys.version.replace("\n", " ")}
    for package in packages[1:]:
        try:
            versions[package] = importlib.metadata.version(package)
        except importlib.metadata.PackageNotFoundError:
            versions[package] = "not installed"
    return versions


def main() -> None:
    args = configure_paths(parse_arguments())
    configure_logging(args.table_dir)

    logging.info("Project root: %s", args.project_root)
    logging.info("Raw loom directory: %s", args.raw_dir)
    logging.info("Metadata: %s", args.metadata)

    metadata = load_metadata(args.metadata)
    if args.validate_files_only:
        validate_input_files(args.raw_dir, metadata)
        logging.info("File-level RNA-velocity input validation completed successfully.")
        return

    adata, matching = load_and_subset_looms(args.raw_dir, metadata)
    logging.info(
        "Combined dataset before velocity preprocessing: %d spots x %d genes",
        adata.n_obs,
        adata.n_vars,
    )

    add_unspliced_qc(adata)
    export_qc(
        adata,
        matching,
        args.table_dir,
        args.diagnostic_unspliced_threshold,
    )

    used_existing_umap = use_existing_umap_if_available(adata)
    preprocess_and_run_velocity(adata, args)

    adata.uns["RNA_velocity_parameters"] = {
        "domains": DOMAIN_ORDER,
        "samples": list(SAMPLE_TO_SUFFIX),
        "min_shared_counts": args.min_shared_counts,
        "n_top_genes": args.n_top_genes,
        "n_pcs": args.n_pcs,
        "n_neighbors": args.n_neighbors,
        "max_iter": args.max_iter,
        "n_jobs": args.n_jobs,
        "velocity_mode": "dynamical",
        "harmony_used": args.use_harmony,
        "precomputed_umap_used": used_existing_umap,
    }

    output_h5ad = (
        args.processed_dir
        / "maize_shoot_SAM_P1_P2_P3_P4_P5_scvelo_dynamical.h5ad"
    )
    adata.write_h5ad(output_h5ad, compression="gzip")
    adata.obs.to_csv(args.table_dir / "RNA_velocity_spot_metadata.csv")
    generate_figures(adata, args.figure_dir)

    with (args.table_dir / "python_package_versions.json").open("w", encoding="utf-8") as handle:
        json.dump(package_versions(), handle, indent=2)

    logging.info("RNA velocity analysis completed.")
    logging.info("Processed AnnData: %s", output_h5ad)
    logging.info("Figures: %s", args.figure_dir)
    logging.info("Tables and log: %s", args.table_dir)


if __name__ == "__main__":
    main()
