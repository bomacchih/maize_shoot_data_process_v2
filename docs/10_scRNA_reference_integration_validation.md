# Validation of `10_scRNA_reference_integration`

Validation date: 2026-08-30

## Outcome

The six R scripts in `scripts/R/10_scRNA_reference_integration/` parse successfully, use project-relative paths, and have consistent input/output handoffs. The available processed objects and required R packages were detected, and the lightweight execution paths for scripts 00–03 completed successfully. Script 04 had already completed Seurat label transfer and section-wise SPOTlight deconvolution; its saved output and restart checkpoint passed the integrity tests below. Script 05 completed against that saved output and regenerated all requested plots.

## Script checks

- `00_prepare_SCINA_marker_table.R`: completed using the project marker-list RDS and produced a validated marker table.
- `01_prepare_maize_scRNA_reference_SCT_Harmony_Seurat_v5.R`: correctly detected the existing final reference and skipped rebuilding it because `overwrite_existing_output` was `FALSE`.
- `02_plot_scRNA_reference_QC_and_Harmony.R`: completed from `sce_ref.rds` and generated panels A and B plus the combined figure.
- `03_scRNA_celltype_annotation_SCINA_Seurat_v5.R`: its plotting/validation path completed from the lightweight reference. The full SCINA result in `sc_merged_filter_SCT2_inte_SCINA.rds` was inspected without rerunning the expensive annotation.
- `04_map_scRNA_to_Visium_SPOTlight_Seurat_v5.R`: the hard-label transfer, 54 section-level deconvolutions, disk checkpoint recovery, metadata transfer, and serialization-safe output path completed. The UMAP reduction selector now recognizes `umap_harmony`.
- `05_plot_SPOTlight_celltypes_Visium_Seurat_v5.R`: completed from the mapped Visium RDS and generated Figure B, Figure C, selected UMAPs, threshold plots, spatial-map PDFs, and a session-information record.

## Saved-object integrity checks

Mapped object:

`data/processed/maize_shoot_14samples_celltype_mapped_SPOTlight_seurat_v5.rds`

- Class: Seurat
- Features: 23,048
- Spots: 20,090
- Assays: `RNA`, `SCT`
- Reductions: `pca`, `umap_pca`, `harmony`, `umap_harmony`
- Spatial images: 14, one for each capture area
- SPOTlight proportion fields: 12
- Non-finite proportion values: 0
- Per-spot proportion sums: exactly 1 for all 20,090 spots
- Active identity levels: 33
- Active identities match `harmony_clusters`: yes
- Checkpoint dimensions: 20,090 spots × 12 reference types
- Checkpoint spot and cell-type names match the mapped object: yes

## Marker and package checks

The prepared SCINA marker table contains 3,997 marker assignments for 14 cell types, with no duplicated gene–cell-type pairs.

Installed versions used for this validation include:

- R 4.6.1
- Seurat 5.5.1
- SeuratObject 5.4.0
- SingleCellExperiment 1.34.0
- SpatialExperiment 1.22.0
- scDblFinder 1.26.7
- glmGamPoi 1.24.0
- harmony 2.0.5
- SCINA 1.2.0
- SPOTlight 1.16.0
- presto 1.0.0
- scuttle 1.22.0
- scran 1.40.0
- scrapper 1.6.3

## Annotation-source warning

The deposited scRNA reference contains two relevant annotation fields, and they are not interchangeable:

- `celltype_scina` is the newly calculated automated SCINA result. It contains 12 non-Unknown classes, including `G2_M_phase` and `S_phase`, but has no cells assigned to `Leaf_guard_cell` or `Pavement_cell_A`.
- `celltype_scina_histo` is the curated 12-tissue-cell-type annotation used for the original Figure B/C representation. It contains `Leaf_guard_cell` and `Pavement_cell_A` and excludes the two cell-cycle labels.

Script 04 now exposes `reference_celltype_column` so this choice is explicit. The currently saved mapped RDS was produced with `celltype_scina`. Reproducing the original 12-tissue-cell-type figure requires setting `reference_celltype_column <- "celltype_scina_histo"` and rebuilding the SPOTlight checkpoint and mapped output under a distinct output name or after intentionally archiving the automated-SCINA result.
