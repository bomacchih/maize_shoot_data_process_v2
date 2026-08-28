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
- one additional processed RDS used as a plotting fallback.

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
