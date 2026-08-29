# 02. Load Space Ranger outputs

This workflow loads one Space Ranger capture-area output from `data/raw/<sample_id>/outs/` into either STUtility or Seurat v5. Enter `sample_id` manually in the selected script and run it from the repository root.

## Scripts

- [`01_load_single_sample_stutility.R`](../scripts/R/02_load_spaceranger/01_load_single_sample_stutility.R)
- [`02_load_single_sample_Seurat_v5.R`](../scripts/R/02_load_spaceranger/02_load_single_sample_Seurat_v5.R)

The STUtility route saves `data/processed/<sample_id>_stutility.rds`. The Seurat v5 route saves `data/processed/<sample_id>_seurat_v5.rds` with `sample_id` and `orig.ident` recorded in the metadata.

## Session information

The script writes the full R session record to:

[`02_load_spaceranger_sessionInfo.txt`](../results/sessionInfo/02_load_spaceranger_sessionInfo.txt)

```r
sessionInfo()
```
