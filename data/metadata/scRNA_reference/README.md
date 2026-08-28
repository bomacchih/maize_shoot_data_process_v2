# scRNA cell-type marker input

`SCINA_marker_table.csv` is the analysis-ready marker file used to annotate the
combined maize scRNA reference before transferring/deconvolving those cell
types in Visium with SPOTlight.

## Source comparison

The curated source, `marker_list2.rds` (MD5
`497fe45631f35884eff358e98bbdd56f`), was compared with **Supplementary Table
9** of `Supplementary_Tables_20251104.xlsx` (MD5
`3d0c20648b193de5f1919ad7127f3b70`). The workbook was treated as data only.

| Check | Supplementary Table 9 | Curated CSV |
|---|---:|---:|
| Raw gene–cell-type rows | 8,080 | 3,997 |
| Repeated gene–cell-type rows | 1,227 | 0 |
| Unique genes | 5,289 | 3,997 |
| Genes assigned to multiple cell types | 1,195 | 0 |
| Cell types | 14 | 14 |

Every curated marker is present in Supplementary Table 9 and is assigned to
only one cell type there. The curated file preserves the Supplementary Table 9
row order after filtering. Supplementary Table 9 contains 4,094 cell-type-
exclusive genes, 97 of which are not in `marker_list2.rds`: 71 are absent from
the combined scRNA feature space, while 26 occur in both the combined scRNA
and representative Visium feature spaces. The 26 were not restored because
the reason for their prior curated exclusion is not recorded; adding them
should be treated as a sensitivity analysis rather than silently changing the
published marker definition.

All 3,997 retained markers occur in both
`sc_merged_filter_SCT2_inte.rds` and the representative Visium object
`UL01_seurat_v5.rds`.

## File schema

| Column | Meaning |
|---|---|
| `gene_id` | Stable maize gene identifier |
| `cell_type` | Analysis-safe cell-type label; punctuation is normalized to underscores |
| `marker_rank` | Position in the curated source list, not an effect-size ranking |
| `source_cell_type` | Original label in `marker_list2.rds` |

The SCINA workflow accepts this file directly. It currently retains at most
100 markers per cell type; `Leaf_rim` has only 9 markers and therefore merits
specific biological review. SPOTlight does not consume this CSV directly: it
uses the SCINA-annotated scRNA identities to calculate reference markers and
map/deconvolve cell types in Visium.

Reproduce the source comparison with
`scripts/python/compare_supplementary_table9_markers.py`. Validate expression-
feature overlap with `tests/R/validate_marker_list_for_mapping.R`.
