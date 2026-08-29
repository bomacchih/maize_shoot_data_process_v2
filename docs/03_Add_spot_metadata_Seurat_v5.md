# 03. Add spot metadata to Seurat v5 objects

These interactive RStudio scripts update a loaded Seurat object without creating a duplicate global object. The metadata CSV must contain `Barcode` in its first column and one or more metadata fields in the remaining columns.

## Scripts

- [`00_fix_SeuratObject_version_issue.R`](../scripts/R/03_add_metadata_seurat5/00_fix_SeuratObject_version_issue.R) updates older Seurat/Visium image objects when required.
- [`01_add_spot_metadata_from_csv_delete_NA.R`](../scripts/R/03_add_metadata_seurat5/01_add_spot_metadata_from_csv_delete_NA.R) imports metadata and removes unmatched spots.
- [`02_add_spot_metadata_from_csv_without_delete_NA.R`](../scripts/R/03_add_metadata_seurat5/02_add_spot_metadata_from_csv_without_delete_NA.R) imports metadata while retaining unmatched spots as `NA`.

The scripts preserve the selected object name and verify that metadata import does not unexpectedly change the active identities.

## Session information

The script writes the full R session record to:

[`03_add_metadata_seurat5_sessionInfo.txt`](../results/sessionInfo/03_add_metadata_seurat5_sessionInfo.txt)

```r
sessionInfo()
```
