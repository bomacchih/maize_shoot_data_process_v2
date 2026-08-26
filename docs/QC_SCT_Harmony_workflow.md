# Quality control, SCTransform, PCA, and Harmony integration

This workflow reproduces the quality-control and integration steps described in the Bio-protocol manuscript for the 14 maize shoot Visium capture areas:

`UL01`, `UL02`, `UL04`, `VR01`-`VR04`, `DQ01`-`DQ04`, and `DQ06`-`DQ08`.

UL03 and DQ05 are not included because individual Seurat RDS files were not provided for these capture areas.

## Input files

Run all scripts from the repository root. The expected inputs are:

```text
data/
├── processed/
│   ├── UL01_seurat_v5.rds
│   ├── ...
│   └── DQ08_seurat_v5.rds
└── reference/
    ├── maize_mitochondrial_genes.txt
    └── maize_plastid_genes.txt
```

The mitochondrial and plastid lists contain 73 and 21 genes, respectively. They were recovered from the original analysis materials, and every listed gene is present in the current maize feature set. Explicit lists are required because these organelle genes use names such as `nad4L-XLOC-043543` and `psbC-XLOC-043729`; a generic `^MT-` or `^CP-` pattern would not identify them.

## 1. Generate QC metrics and plots

Run:

```r
source("scripts/R/04_merge_integration/00_QC_14_samples_Seurat_v5.R")
```

The script calculates the following metrics for every tissue-covered spot:

- `nFeature_RNA`: number of detected genes.
- `nCount_RNA`: total UMI count.
- `percent.mito`: percentage of UMIs assigned to the mitochondrial gene list.
- `percent.pltd`: percentage of UMIs assigned to the plastid gene list.

It also summarizes each gene using:

- `total_umi`: total UMI count across all 14 datasets.
- `detected_spots`: number of tissue spots in which the gene is detected.
- `pass_minimum_100_reads`: whether the gene has at least 100 total reads.

The manuscript reports that mitochondrial reads were below 4% and plastid reads were below 1% in most tissue-covered spots. These values are descriptive observations from the original dataset, not universal filtering thresholds. Inspect the distributions for each new dataset before deciding whether spot-level filtering is necessary.

The script produces:

```text
results/
├── figures/
│   ├── QC_spot_gene_histograms.png
│   ├── QC_violin_by_sample.png
│   └── QC_features_vs_counts.png
└── tables/
    ├── QC_sample_summary.csv
    └── QC_gene_summary.csv
```

`QC_spot_gene_histograms.png` contains:

- A: unique genes per tissue spot.
- B: total UMI counts per tissue spot.
- C: total UMI counts per gene on a log10 scale.
- D: total tissue spots per gene.

`QC_violin_by_sample.png` shows the four spot-level QC metrics separately for all 14 capture areas. `QC_features_vs_counts.png` checks the relationship between detected genes and total UMI counts.

### QC distributions across spots and genes

![QC distributions across tissue spots and genes](../results/figures/QC_spot_gene_histograms.png)

**Figure 1. Quality-control distributions across tissue spots and genes.** (A) Number of unique genes detected per tissue spot. (B) Total mapped UMI counts per tissue spot. (C) Total mapped UMI counts per gene on a log10 scale. (D) Number of tissue spots in which each gene was detected.

### QC metrics for individual capture areas

![Violin plots of QC metrics for the 14 capture areas](../results/figures/QC_violin_by_sample.png)

**Figure 2. Spot-level QC metrics across the 14 capture areas.** Violin plots show the distributions of detected genes (`nFeature_RNA`), total UMI counts (`nCount_RNA`), mitochondrial-read percentages (`percent.mito`), and plastid-read percentages (`percent.pltd`) for each sample.

### Relationship between detected genes and UMI counts

![Relationship between detected genes and total UMI counts](../results/figures/QC_features_vs_counts.png)

**Figure 3. Relationship between library complexity and detected features.** Each point represents one tissue-covered spot. The plot is used to identify spots with unusually low gene detection, low UMI counts, or atypical count-to-feature relationships.

## 2. Apply the gene-level QC rule and integrate datasets

Run:

```r
source("scripts/R/04_merge_integration/01_merge_14_samples_SCT_Harmony_Seurat_v5.R")
```

The integration script performs the following steps:

1. Loads and updates the 14 Seurat v5 objects.
2. Standardizes the expression assay name to `RNA` when necessary.
3. Adds `sample_id`, `orig.ident`, `percent.mito`, and `percent.pltd` metadata.
4. Merges the objects while retaining the 14 sample-specific RNA count layers.
5. Removes genes with fewer than 100 total reads across the combined dataset, as specified in the Bio-protocol.
6. Normalizes each sample layer using SCTransform v2 with 3,000 variable features.
7. Calculates 50 PCs for the scree/elbow diagnostic.
8. Records the automatic elbow candidate but retains PCs 1-30 for the Bio-protocol analysis.
9. Recalculates the final PCA with 30 dimensions.
10. Generates an uncorrected PCA-based UMAP.
11. Runs Seurat v5 `IntegrateLayers()` with `HarmonyIntegration` using PCs 1-30.
12. Generates a Harmony-based UMAP and a before-versus-after integration figure.
13. Validates and saves the combined Seurat object.

## PCA selection

The elbow diagnostic evaluates PCs 1-50. The plot contains two reference lines:

- Blue dotted line: automatically detected strongest elbow.
- Red dashed line: 30 PCs retained for the Bio-protocol workflow.

Thirty PCs are retained to match the original analysis and provide a common dimensional space for integration, UMAP, and later clustering. The automatic elbow is reported separately for transparency. When applying this workflow to another dataset, users should inspect the scree plot and evaluate whether a different number of PCs better preserves biological structure without adding low-variance noise.

PCA outputs are:

```text
results/
├── figures/PCA_elbow_plot.png
├── tables/PCA_variance_explained.csv
└── logs/PCA_selection.txt
```

![PCA elbow plot used to select the retained dimensions](../results/figures/PCA_elbow_plot.png)

**Figure 4. PCA elbow diagnostic.** Variance explained is shown for PCs 1-50. The blue dotted line marks the automatically detected strongest elbow, whereas the red dashed line marks PC 30, the number of PCs retained to reproduce the Bio-protocol workflow.

## Before and after Harmony integration

`results/figures/PCA_Harmony_before_after_UMAP.png` compares:

- A: UMAP constructed from the uncorrected PCA space.
- B: UMAP constructed from the SCTransform-normalized, Harmony-corrected space.

Both panels are colored by `sample_id`. Better mixing after Harmony is consistent with reduced sample-associated variation, but biological and anatomical structure should also be checked before accepting the integration.

![UMAP before and after SCTransform and Harmony integration](../results/figures/PCA_Harmony_before_after_UMAP.png)

**Figure 5. Sample distributions before and after Harmony integration.** (A) UMAP generated from the uncorrected PCA representation after SCTransform normalization. (B) UMAP generated from the Harmony-integrated representation using PCs 1-30. Spots are colored by capture-area identity. Increased mixing of samples after integration is consistent with reduced sample-associated effects.

## Final output

The integrated object is saved as:

```text
data/processed/maize_shoot_14samples_SCT_harmony_seurat_v5.rds
```

The expected final reductions are:

- `pca`: 30 dimensions.
- `harmony`: 30 dimensions.
- `umap_pca`: UMAP before integration.
- `umap_harmony`: UMAP after integration.

The object retains `RNA` and `SCT` assays, the 14 sample identities, and the organelle QC metadata. The obsolete `nCount_Spatial` and `nFeature_Spatial` columns are not used.

### Verified result for the current 14-sample test run

- Tissue spots: 23,160.
- Raw merged feature set: 40,109 genes.
- Genes retained after the 100-read rule: 23,550.
- Automatic strongest elbow candidate: PC 6.
- PCs retained for the Bio-protocol analysis: PCs 1-30.
- PCA dimensions in the saved object: 30.
- Harmony dimensions in the saved object: 30.
- Median mitochondrial percentage: 0.0185%; maximum: 6.61%.
- Median plastid percentage: 0%; maximum: 0.91%.
- Final Seurat object validation: passed.

## Notes for downstream analyses

- Tissue spots that overlap two or more anatomical domains should be excluded from domain-level comparisons, as described in the manuscript. They do not need to be removed from the full object before integration.
- Do not interpret spots from the same section or capture area as independent biological replicates.
- For pseudobulk or differential analyses, aggregate and model counts at the biological-replicate level rather than treating individual spots as replicates.
- Use the Harmony reduction for integrated UMAP, neighbor finding, and clustering. Use raw RNA counts for count-based differential-expression or pseudobulk models.

---

## Session information

Run this command at the end of the analysis to record the R version, operating system, and loaded package versions:

```r
sessionInfo()
```
