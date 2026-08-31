# Published processed data

Download the public files from:

- Zenodo record: <https://zenodo.org/records/22058284>
- DOI: <https://doi.org/10.5281/zenodo.22058284>

Place the downloaded files directly in this directory, keeping their deposited
filenames. The record was checked through the Zenodo API on 2026-08-28 and
contains:

- 14 Space Ranger `web_summary.html` files;
- 14 Loupe Browser `.cloupe` files;
- 14 Velocyto `.loom` files;
- 14 per-sample Seurat `.rds` objects;
- the integrated Visium object
  `maize_shoot_14samples_SCT_harmony_seurat_v5.rds`;
- the integrated scRNA reference `sc_merged_filter_SCT2_inte.rds`; and
- the embryonic-leaf object
  `XGE202122_S5_subset_embleaf_harmony_join.rds`, which defines the 6,392-spot
  SAM–P5 analysis scope and provides the authoritative `umap.harmony`
  coordinates for both the Monocle and SPOTlight workflows.

These objects use different historical metadata conventions:

| Seurat object | Biological replicate | Physical section |
|---|---|---|
| `maize_shoot_14samples_SCT_harmony_seurat_v5.rds` | `sample_id` | `section_id` |
| `XGE202122_S5_subset_embleaf_harmony_join.rds` | `sample` | `sample_id` |

Downstream code must resolve the meaning by object rather than transferring
the column name blindly. The RNA-velocity export adds explicit
`biological_replicate` and `section_id` fields without changing either RDS.

The derived file
`XGE202122_S5_subset_embleaf_celltype_mapped_SPOTlight_seurat_v5.rds` is
generated locally by recalculating Seurat label transfer and SPOTlight
deconvolution on that subset; it is intentionally not tracked in Git.

## Important distinction

The web summaries and `.cloupe` files are Space Ranger deliverables, but the
Zenodo record does not include the complete Space Ranger `outs` directories.
In particular, it does not contain `filtered_feature_bc_matrix.h5` or the
spatial image, position, and scalefactor files required by the STUtility and
Seurat loading scripts. Those full inputs remain necessary only when rebuilding
the per-sample Seurat objects from Space Ranger results. Most downstream
workflows can start from the deposited `.rds` objects instead.

The `.loom` files deposited here are the inputs for the RNA-velocity script.
The combined scRNA reference RDS supports downstream SCINA annotation; the
original SRR11943512 and SRR11943513 10x matrices are necessary only to rebuild
that RDS from raw counts.
