# 03. Add spot metadata to Seurat v5 objects

These interactive RStudio scripts update a loaded Seurat object without creating a duplicate global object. The metadata CSV must contain `Barcode`; matching is performed against the original 10x barcode extracted from each Seurat spot name.

## Scripts

- [`00_fix_SeuratObject_version_issue.R`](../scripts/R/03_add_metadata_seurat5/00_fix_SeuratObject_version_issue.R) updates older Seurat/Visium image objects when required.
- [`01_add_spot_metadata_from_csv_delete_NA.R`](../scripts/R/03_add_metadata_seurat5/01_add_spot_metadata_from_csv_delete_NA.R) requires `orig.ident`, `section_id`, `domains`, `sample_id`, `section`, `cca_clusters`, `seurat_clusters`, `sample_domain`, `ms_ve`, and `domain_section`; it imports these fields and removes unmatched spots.
- [`02_add_spot_metadata_from_csv_without_delete_NA.R`](../scripts/R/03_add_metadata_seurat5/02_add_spot_metadata_from_csv_without_delete_NA.R) imports every CSV column except `Barcode`. Unmatched spots are retained with an existing value when present, or `NA` for a newly added field.

Both scripts preserve the selected object name and verify that metadata import does not unexpectedly change the active identities. For the filtering script, identities are compared only among the retained matched spots.

## Session information

The script writes the full R session record to:

[`03_add_metadata_seurat5_sessionInfo.txt`](../results/sessionInfo/03_add_metadata_seurat5_sessionInfo.txt)

```r
sessionInfo()
```
