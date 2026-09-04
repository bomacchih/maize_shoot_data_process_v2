# 10. scRNA reference and SAM–P5 Visium mapping validation

Validation date: 2026-09-01

This report records the current subset-based scRNA-to-Visium workflow. The
canonical procedure and figures are documented in
[`10_scRNA_to_Visium_SPOTlight_mapping.md`](10_scRNA_to_Visium_SPOTlight_mapping.md).
Older full-object checkpoints and tables may remain in `results/`, but they are
legacy outputs and are not the inputs used by the current scripts.

## Current analysis scope

- Annotated scRNA reference:
  `data/processed/sc_merged_filter_SCT2_inte_SCINA.rds`
- Visium query and coordinate source:
  `data/processed/XGE202122_S5_subset_embleaf_harmony_join.rds`
- Query membership: 6,392 spots assigned to `SAM`, `P1_P2`, `P3`, `P4`, or
  `P5`
- UMAP coordinates: the exact `umap.harmony` embedding stored in the deposited
  SAM–P5 query object
- Current mapped output:
  `data/processed/XGE202122_S5_subset_embleaf_celltype_mapped_SPOTlight_seurat_v5.rds`

Both Seurat anchor transfer and SPOTlight deconvolution are recalculated on
this same 6,392-spot subset. The script does not reuse hard labels calculated
on the complete 20,090-spot object unless that behavior is explicitly enabled.

## Current completed outputs

The subset deconvolution table contains all 6,392 spots. The section diagnostic
records nine section groups (`Section0`–`Section8`), 28,750 shared genes, 3,000
HVGs, and 1,200 retained marker entries per completed section. Per-spot
SPOTlight proportions are normalized to sum to 1.

The current Seurat-versus-SPOTlight agreement summary contains 6,392 evaluable
spots. Agreement is a diagnostic comparison between a hard maximum-posterior
label and a soft mixture model; it is not an accuracy estimate or a filtering
criterion.

Important current files are:

- `results/tables/10_scRNA_Visium_mapping/SPOTlight_subset_embleaf_spot_celltype_proportions.csv`
- `results/tables/10_scRNA_Visium_mapping/SPOTlight_subset_embleaf_section_diagnostics.csv`
- `results/tables/10_scRNA_Visium_mapping/SPOTlight_subset_embleaf_vs_Seurat_agreement_summary.csv`
- `results/figures/10_scRNA_Visium_mapping/selected_celltypes/Figure_B_VR03_section_scatterpies.png`
- `results/figures/10_scRNA_Visium_mapping/selected_celltypes/Figure_C_SPOTlight_proportion_UMAPs.png`

## Annotation-field choice

The reference may contain both `celltype_scina` and `celltype_scina_histo`.
Script 04 exposes `reference_celltype_column` so this selection can be made
explicitly. Use `celltype_scina` for the newly calculated automated SCINA
annotation. Use `celltype_scina_histo` only when intentionally reproducing the
curated 12-tissue-cell-type Figure B/C representation. Do not describe an
output as one annotation scheme unless that field was selected for the run.

## Session information

The mapping and plotting scripts write their full R session records to:

- [`12_scRNA_Visium_mapping_sessionInfo.txt`](../results/sessionInfo/12_scRNA_Visium_mapping_sessionInfo.txt)
- [`12_scRNA_Visium_plotting_sessionInfo.txt`](../results/sessionInfo/12_scRNA_Visium_plotting_sessionInfo.txt)

```r
sessionInfo()
```
