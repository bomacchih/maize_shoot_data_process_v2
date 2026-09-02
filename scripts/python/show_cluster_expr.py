#!/usr/bin/env python3
"""Plot a TO-GCN transcription-factor expression heatmap.

The command accepts an expression matrix, a TF-to-level table, a PNG output,
and a CSV output.  Source level numbers can be shifted for display without
modifying the input data.
"""

from __future__ import annotations

import argparse
from collections.abc import Sequence
from pathlib import Path
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import leaves_list, linkage
from scipy.spatial.distance import pdist


def build_parser() -> argparse.ArgumentParser:
    """Build and return the command-line parser."""
    parser = argparse.ArgumentParser(
        description=(
            "Order transcription factors within TO-GCN levels, calculate "
            "gene-wise z-scores, and write a heatmap plus scaled-expression CSV."
        )
    )
    parser.add_argument(
        "expression_tsv",
        type=Path,
        help="Tab-separated expression matrix with gene IDs in the first column.",
    )
    parser.add_argument(
        "level_csv",
        type=Path,
        help="CSV assigning TF gene IDs (first column) to levels (last column).",
    )
    parser.add_argument("output_png", type=Path, help="Output heatmap PNG path.")
    parser.add_argument(
        "output_csv",
        type=Path,
        help="Output gene-wise z-scored expression CSV path.",
    )
    parser.add_argument(
        "--level-offset",
        type=int,
        default=1,
        help=(
            "Integer subtracted from source levels when labels are displayed "
            "(default: 1, so source levels 2-14 become L1-L13)."
        ),
    )
    return parser


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    """Parse command-line arguments from *argv* or the process command line."""
    return build_parser().parse_args(argv)


def display_level(source_level: int, level_offset: int) -> str:
    """Return the plot label for a source level after applying the offset."""
    return f"L{source_level - level_offset}"


def _validated_gene_ids(values: pd.Index | pd.Series, table_name: str) -> pd.Index:
    """Return stripped gene IDs after validating blanks and duplicates."""
    series = pd.Series(values, dtype="object")
    if series.isna().any():
        raise ValueError(f"{table_name} contains blank gene IDs.")
    gene_ids = pd.Index(series.astype(str).str.strip())
    if (gene_ids == "").any():
        raise ValueError(f"{table_name} contains blank gene IDs.")
    if gene_ids.has_duplicates:
        raise ValueError(f"{table_name} contains duplicate gene IDs.")
    return gene_ids


def load_inputs(
    expression_path: Path,
    levels_path: Path,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Load, validate, and intersect expression and TF-level input tables.

    The returned expression table preserves expression-file row order. The
    returned level table is indexed by TF gene ID in level-file row order and
    has one integer column named ``source_level``.
    """
    expression = pd.read_csv(expression_path, sep="\t", index_col=0)
    if expression.shape[1] == 0:
        raise ValueError("Expression table must contain at least one data column.")
    expression.index = _validated_gene_ids(expression.index, "Expression table")

    numeric_expression = expression.apply(pd.to_numeric, errors="coerce")
    if numeric_expression.isna().any().any():
        raise ValueError("Expression table contains nonnumeric expression values.")
    if not np.isfinite(numeric_expression.to_numpy(dtype=float)).all():
        raise ValueError("Expression table contains non-finite expression values.")
    expression = numeric_expression.astype(float)

    raw_levels = pd.read_csv(levels_path)
    if raw_levels.shape[1] < 2:
        raise ValueError("TF-level table must contain gene-ID and level columns.")
    gene_ids = _validated_gene_ids(raw_levels.iloc[:, 0], "TF-level table")
    numeric_levels = pd.to_numeric(raw_levels.iloc[:, -1], errors="coerce")
    if numeric_levels.isna().any() or not np.isfinite(
        numeric_levels.to_numpy(dtype=float)
    ).all():
        raise ValueError("TF-level table must contain numeric integer levels.")
    if not np.equal(numeric_levels, np.floor(numeric_levels)).all():
        raise ValueError("TF-level table must contain numeric integer levels.")

    levels = pd.DataFrame(
        {"source_level": numeric_levels.to_numpy(dtype=int)},
        index=gene_ids,
    )
    expression_gene_count = len(expression)
    level_gene_count = len(levels)
    common_gene_ids = expression.index.intersection(levels.index, sort=False)
    if len(common_gene_ids) == 0:
        raise ValueError("No gene IDs are shared by the two input tables.")

    expression = expression.loc[common_gene_ids]
    levels = levels.loc[levels.index.intersection(expression.index, sort=False)]
    expression.attrs["input_gene_count"] = expression_gene_count
    levels.attrs["input_gene_count"] = level_gene_count
    return expression, levels


def order_genes_by_level(
    expression: pd.DataFrame,
    levels: pd.DataFrame,
) -> tuple[list[str], list[tuple[int, list[str]]]]:
    """Order genes within ascending numeric levels by correlation clustering.

    Constant rows have undefined Pearson correlation distance. They are sorted
    by gene ID after the variable genes in their level. A level containing one
    variable gene requires no clustering.
    """
    ordered_genes: list[str] = []
    groups: list[tuple[int, list[str]]] = []

    for source_level in sorted(levels["source_level"].unique()):
        level_gene_ids = levels.index[
            levels["source_level"] == source_level
        ].tolist()
        level_values = expression.loc[level_gene_ids].to_numpy(dtype=float)
        constant_mask = np.ptp(level_values, axis=1) == 0

        variable_gene_ids = [
            gene_id
            for gene_id, is_constant in zip(
                level_gene_ids, constant_mask, strict=True
            )
            if not is_constant
        ]
        constant_gene_ids = sorted(
            gene_id
            for gene_id, is_constant in zip(
                level_gene_ids, constant_mask, strict=True
            )
            if is_constant
        )

        if len(variable_gene_ids) > 1:
            variable_values = expression.loc[variable_gene_ids].to_numpy(
                dtype=float
            )
            distances = np.clip(
                pdist(variable_values, metric="correlation"), 0.0, 2.0
            )
            if not np.isfinite(distances).all():
                raise ValueError(
                    f"Could not calculate finite correlation distances for "
                    f"source level {source_level}."
                )
            tree = linkage(
                distances,
                method="average",
                optimal_ordering=True,
            )
            variable_gene_ids = [
                variable_gene_ids[position] for position in leaves_list(tree)
            ]

        level_order = variable_gene_ids + constant_gene_ids
        ordered_genes.extend(level_order)
        groups.append((int(source_level), level_order))

    return ordered_genes, groups


def row_zscore(expression: pd.DataFrame) -> pd.DataFrame:
    """Return population z-scores across columns for each expression row.

    Rows with zero variance are returned as zeroes because they have no
    within-gene change to display.
    """
    values = expression.to_numpy(dtype=float)
    means = values.mean(axis=1, keepdims=True)
    standard_deviations = values.std(axis=1, ddof=0, keepdims=True)
    scaled = np.divide(
        values - means,
        standard_deviations,
        out=np.zeros_like(values, dtype=float),
        where=standard_deviations != 0,
    )
    return pd.DataFrame(scaled, index=expression.index, columns=expression.columns)


def plot_heatmap(
    scaled_expression: pd.DataFrame,
    groups: list[tuple[int, list[str]]],
    level_offset: int,
    output_path: Path,
) -> None:
    """Write a level-annotated heatmap of transposed scaled expression."""
    figure, axis = plt.subplots(figsize=(16, 4.5))
    image = axis.imshow(
        scaled_expression.to_numpy(dtype=float).T,
        aspect="auto",
        interpolation="nearest",
        cmap="seismic",
        vmin=-1.5,
        vmax=1.5,
    )

    tick_positions: list[float] = []
    tick_labels: list[str] = []
    group_start = 0
    for group_number, (source_level, gene_ids) in enumerate(groups):
        group_stop = group_start + len(gene_ids)
        tick_positions.append((group_start + group_stop - 1) / 2)
        tick_labels.append(display_level(source_level, level_offset))
        if group_number < len(groups) - 1:
            axis.axvline(group_stop - 0.5, color="gray", linewidth=0.7)
        group_start = group_stop

    axis.set_xticks(tick_positions, labels=tick_labels, fontsize=12)
    axis.set_yticks(
        np.arange(len(scaled_expression.columns)),
        labels=scaled_expression.columns,
        fontsize=12,
    )
    axis.set_xlabel("TO-GCN level", fontsize=12)
    axis.set_ylabel("Developmental domain", fontsize=12)
    axis.tick_params(axis="x", length=0, pad=8)

    colorbar = figure.colorbar(image, ax=axis, pad=0.015)
    colorbar.set_label("Z-score", fontsize=13)
    colorbar.ax.tick_params(labelsize=11)

    figure.tight_layout()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(output_path, dpi=500, bbox_inches="tight")
    plt.close(figure)


def main(argv: Sequence[str] | None = None) -> int:
    """Run the command-line program."""
    args = parse_args(argv)
    try:
        expression, levels = load_inputs(args.expression_tsv, args.level_csv)
        ordered_gene_ids, groups = order_genes_by_level(expression, levels)
        ordered_expression = expression.loc[ordered_gene_ids]
        scaled_expression = row_zscore(ordered_expression)

        args.output_csv.parent.mkdir(parents=True, exist_ok=True)
        scaled_expression.to_csv(args.output_csv)
        plot_heatmap(
            scaled_expression,
            groups,
            args.level_offset,
            args.output_png,
        )
    except (OSError, ValueError, pd.errors.ParserError) as error:
        print(f"error: {error}", file=sys.stderr)
        return 1

    print(f"Matched genes: {len(ordered_gene_ids):,}")
    print(
        "Expression-only genes excluded: "
        f"{expression.attrs['input_gene_count'] - len(ordered_gene_ids):,}"
    )
    print(
        "Level-only genes excluded: "
        f"{levels.attrs['input_gene_count'] - len(ordered_gene_ids):,}"
    )
    print(f"Wrote heatmap: {args.output_png}")
    print(f"Wrote scaled expression: {args.output_csv}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
