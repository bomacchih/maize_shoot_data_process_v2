# Processed datasets to download

Large processed objects are intentionally excluded from Git. Download the
required RDS, loom, `.cloupe`, and Space Ranger web-summary files from
[Zenodo record 22058284](https://zenodo.org/records/22058284), then place them
in the pipeline locations documented in [`data/DATASETS.md`](../DATASETS.md).

- Place Seurat RDS objects, web summaries, and `.cloupe` files in
  `data/processed/`.
- Place the 14 Velocyto loom files in `data/raw/loom/`.
- Do not rename files unless the corresponding script configuration is also
  updated.

This directory is only a Git-tracked download reminder; analysis scripts do
not read data from it.
