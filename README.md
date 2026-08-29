# Maize shoot spatial transcriptomics analysis

Reproducible workflows for processing, integrating, and analyzing spatial transcriptomes from developing maize shoots. This repository contains the revised R, Python, and shell pipelines associated with the Bio-protocol manuscript and the original *Plant Biotechnology Journal* study.

The workflow covers 10x Genomics Visium processing, Seurat v5 quality control and integration, anatomical annotation, pseudobulk analysis, developmental expression trends, RNA velocity, Monocle 3 pseudotime, maize scRNA-seq reference construction, SCINA annotation, and SPOTlight deconvolution.

> **Project status:** the repository is being prepared for manuscript revision and reproducible release. Scripts use project-relative paths, but some raw-data workflows still require external files or local configuration. Check [data/DATASETS.md](data/DATASETS.md) before running an analysis.

## Study design

The principal Visium dataset contains 14 capture areas in the following fixed order:

```text
UL01, UL02, UL04,
VR01, VR02, VR03, VR04,
DQ01, DQ02, DQ03, DQ04, DQ06, DQ07, DQ08
```

Curated analyses retain tissue-covered, non-overlapping spots assigned to one of seven structural domains:

```text
SAM → P1_P2 → P3 → P4 → P5 / (co_v + coleoptile)
```

The first five domains are used for developmental trajectory, RNA-velocity, and pseudotime analyses. The sample order also determines merged barcode suffixes (`_1_1` through `_1_14`); changing it will break correspondence with the curated metadata.

## Workflow overview

| Step | Analysis | Main implementation |
|---:|---|---|
| 1 | Build the maize reference and run Space Ranger | [`scripts/shell/01_spaceranger_processing/`](scripts/shell/01_spaceranger_processing/) |
| 2 | Load one capture area into STUtility or Seurat v5 | [`scripts/R/02_load_spaceranger/`](scripts/R/02_load_spaceranger/) |
| 3 | Add manually curated spot metadata | [`scripts/R/03_add_metadata_seurat5/`](scripts/R/03_add_metadata_seurat5/) |
| 4 | QC, merge 14 samples, run SCTransform, PCA, Harmony, UMAP, and clustering | [`scripts/R/04_merge_integration/`](scripts/R/04_merge_integration/) |
| 5 | Find cluster markers and assign 12 tissue supergroups | [`scripts/R/05_marker_analysis/`](scripts/R/05_marker_analysis/) |
| 6 | Generate replicate-level pseudobulk profiles and Figure 13 | [`scripts/R/06_pseudobulk_analysis/`](scripts/R/06_pseudobulk_analysis/) |
| 7 | Cluster developmental expression trends, analyze GO enrichment, and generate Figure 9 | [`scripts/R/07_developmental_trends_GO/`](scripts/R/07_developmental_trends_GO/) |
| 8 | Estimate dynamical RNA velocity with scVelo | [`scripts/python/08_RNA_velocity/`](scripts/python/08_RNA_velocity/) |
| 9 | Estimate pseudotime and trajectories with Monocle 3 | [`scripts/R/09_monocle3_pseudotime/`](scripts/R/09_monocle3_pseudotime/) |
| 10 | Build and annotate the scRNA-seq reference; map cell types to Visium with Seurat and SPOTlight | [`scripts/R/10_scRNA_reference_integration/`](scripts/R/10_scRNA_reference_integration/) |

The numbered directories indicate the intended execution order. Within a directory, run scripts in filename order.

## Repository layout

```text
maize_shoot_data_process_v2/
├── config/                         # Pipeline configuration files
├── data/
│   ├── DATASETS.md                 # Human-readable input and readiness inventory
│   ├── dataset_inventory.csv       # Machine-readable dataset inventory
│   ├── metadata/                   # Curated spot and section metadata
│   ├── processed/                  # Large local outputs; excluded from Git
│   ├── processed_need_to_download/ # Placeholder and download instructions
│   ├── raw/                        # Raw/local inputs; large files excluded from Git
│   └── reference/                  # Small gene and GO reference tables
├── docs/                           # Workflow reports and Markdown sources
├── environment/                    # Environment specifications
├── results/                        # Generated figures, tables, and run information
├── scripts/
│   ├── R/                          # Seurat, Monocle 3, SCINA, and SPOTlight
│   ├── python/                     # scVelo dynamical RNA velocity
│   └── shell/                      # Space Ranger processing
├── tests/                          # Validation utilities
└── workflow/                       # Workflow-level orchestration files
```

Large sequencing files, Space Ranger outputs, loom files, Seurat objects, and other derived binaries are intentionally excluded from normal Git tracking.

## Data availability

| Dataset | Accession |
|---|---|
| Raw Visium sequencing reads | NCBI BioProjects [PRJNA805024](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA805024) and [PRJNA804974](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA804974) |
| Processed spatial transcriptomics data | GEO [GSE196882](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE196882) |
| Maize shoot scRNA-seq reference | SRA [SRR11943512](https://www.ncbi.nlm.nih.gov/sra/?term=SRR11943512) and [SRR11943513](https://www.ncbi.nlm.nih.gov/sra/?term=SRR11943513) |
| Napari interactive-viewer datasets | Zenodo [16933147](https://zenodo.org/records/16933147) |
| RDS, loom, `.cloupe`, and web-summary datasets | Zenodo [22058284](https://zenodo.org/records/22058284) |

See [data/DATASETS.md](data/DATASETS.md) for required filenames, expected locations, availability, external-data locations, and missing inputs. Reviewer links and private repository tokens must never be committed to Git.

## Software requirements

The primary analysis environment uses:

- R 4.x with Seurat v5-compatible packages
- Python 3.10–3.12 for scVelo
- 10x Genomics Space Ranger for reference construction and Visium processing
- Velocyto.py 0.17 series for generating spliced/unspliced loom files

Detailed requirements are recorded in:

- [data/SOFTWARE_REQUIREMENTS.md](data/SOFTWARE_REQUIREMENTS.md) — package roles and audited versions
- [data/R_packages.csv](data/R_packages.csv) — machine-readable R package inventory
- [data/requirements-python.txt](data/requirements-python.txt) — Python dependencies

Install Python dependencies in an isolated environment, for example:

```bash
python -m venv .venv
# Windows: .venv\Scripts\activate
# Linux/macOS: source .venv/bin/activate
python -m pip install -r data/requirements-python.txt
```

Bioconductor and GitHub packages should be installed using the sources documented in [data/SOFTWARE_REQUIREMENTS.md](data/SOFTWARE_REQUIREMENTS.md). Record the final R and Python environments before release.

## Quick start from processed Seurat objects

This is the shortest entry point when the 14 individual Seurat objects have already been downloaded.

1. Clone the repository and enter its root directory.

   ```bash
   git clone https://github.com/bomacchih/maize_shoot_data_process_v2.git
   cd maize_shoot_data_process_v2
   ```

2. Place the 14 files in `data/processed/` using names such as `UL01_seurat_v5.rds` and `DQ08_seurat_v5.rds`.

3. Confirm that `data/metadata/metadata.csv` is present. This table defines the accepted non-overlapping shoot spots and supplies structural-domain annotations.

4. Start R from the repository root and run QC followed by integration:

   ```r
   source("scripts/R/04_merge_integration/00_QC_14_samples_Seurat_v5.R")
   source("scripts/R/04_merge_integration/01_merge_14_samples_SCT_Harmony_Seurat_v5.R")
   ```

5. Inspect the PCA elbow plot and QC outputs. The reference workflow retains 30 PCs because the contribution of later components decreased substantially, but users should adjust this value for their own data.

6. Continue with the numbered downstream scripts required for the desired analysis.

The integration workflow writes the main combined object to:

```text
data/processed/maize_shoot_14samples_SCT_harmony_seurat_v5.rds
```

The object contains two deliberately separate cluster fields:

- `harmony_clusters_recomputed` is calculated from the current Harmony reduction using `FindNeighbors()` and `FindClusters()` and is the appropriate starting point for a new analysis.
- `harmony_clusters` is imported from `metadata.csv` and retained as the active identity to reproduce the published marker analysis and Figure 12 mapping exactly.

Numeric cluster labels are analysis-specific. Users must characterize their recomputed clusters and create their own cluster-to-supergroup mapping rather than transferring the published mapping.

## Starting from Space Ranger outputs

Each capture area should follow this layout:

```text
data/raw/<sample_id>/outs/
├── filtered_feature_bc_matrix.h5
└── spatial/
    ├── tissue_hires_image.png
    ├── tissue_positions.csv
    └── scalefactors_json.json
```

Use one of the loaders in [`scripts/R/02_load_spaceranger/`](scripts/R/02_load_spaceranger/) and manually set the sample ID and input directory. The STUtility loader may require the legacy filename `tissue_positions_list.csv`; consult [data/DATASETS.md](data/DATASETS.md) before loading newer Space Ranger outputs.

The shell scripts contain installation- and project-specific settings for FASTQs, histology images, transcriptome references, slide identifiers, and capture areas. Review those variables before submitting Space Ranger jobs.

## RNA velocity

Place one Velocyto loom file for each of the 14 samples in `data/processed/`. Filenames must contain their sample IDs. Validate filenames and metadata correspondence without running the numerical analysis:

```bash
python scripts/python/08_RNA_velocity/01_scvelo_dynamical_RNA_velocity.py --validate-files-only
```

Run the complete dynamical model with:

```bash
python scripts/python/08_RNA_velocity/01_scvelo_dynamical_RNA_velocity.py
```

The workflow retains SAM, P1_P2, P3, P4, and P5 spots, reports low-unspliced-fraction groups, fits the scVelo dynamical model, and exports an H5AD object, figures, tables, a run log, and Python package versions. Harmony correction is optional through `--use-harmony` and is off by default because it was not specified for RNA velocity in the manuscript method.

## Outputs and reports

Scripts write reproducible products to consistent locations:

```text
data/processed/   RDS, H5AD, and other reusable processed objects
results/figures/  PNG and PDF figures
results/tables/   CSV and compressed matrix exports
results/logs/ and results/sessionInfo/  Run summaries and environment records
docs/             Human-readable Markdown and rendered HTML reports
```

Available reports include:

- [QC, SCTransform, Harmony integration, and clustering](docs/QC_SCT_Harmony_workflow.md)
- [Seurat v5 maize spatial dataset](docs/Maize_data_Seurat_v5.md)
- [Tissue supergroups and Figure 12](docs/Tissue_supergroups_Figure_12_Seurat_v5.md)
- [Structural-domain pseudobulk analysis and Figure 13](docs/Pseudobulk_structural_domains_Figure_13.md)
- [Developmental trends and GO enrichment](docs/Developmental_expression_trends_GO_Figure_9.md)
- [scRNA-seq reference QC and Harmony](docs/scRNA_reference_QC_and_Harmony.md)
- [SCINA cell-type annotation](docs/scRNA_reference_SCINA_celltype_annotation.md)
- [scRNA-to-Visium SPOTlight mapping](docs/scRNA_to_Visium_SPOTlight_mapping.md)

## Reproducibility notes

- Run scripts from the repository root unless a script explicitly documents another entry point.
- Do not reorder the 14 sample IDs without rebuilding barcode-to-metadata mappings.
- Do not treat spots from the same section or capture area as independent biological replicates. The pseudobulk workflow aggregates raw counts by biological replicate and anatomical domain.
- Inspect QC summaries, the PCA elbow plot, batch mixing, marker expression, and anatomical correspondence before accepting downstream results.
- Keep raw counts for pseudobulk and differential analyses; the SCT assay is used for visualization and selected expression-pattern analyses.
- Scripts do not overwrite major processed objects silently where an overwrite guard is implemented.
- Each completed analysis should retain R `sessionInfo()` or the Python package-version record generated by the workflow.

Before changing the repository visibility or preparing a release, run the
[public-release audit and checklist](docs/PUBLIC_RELEASE_CHECKLIST.md).

## Citation

This protocol was used in:

Wu et al. *Plant Biotechnology Journal* (2026). [https://doi.org/10.1111/pbi.70515](https://doi.org/10.1111/pbi.70515)

When using this repository, cite both the Bio-protocol article after publication and the original study. Also cite the primary software used in the relevant analysis, including Seurat, Harmony, STUtility, scVelo, Monocle 3, SCINA, and SPOTlight.

## Questions and issue reports

When reporting a reproducibility problem, include the script name, command or R call, complete error message, relevant input-file layout, and `sessionInfo()` or Python package versions. Do not include reviewer URLs, secure tokens, credentials, or unpublished controlled-access data in a public issue.
