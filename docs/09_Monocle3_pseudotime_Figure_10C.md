# 09. Monocle 3 pseudotime analysis of maize embryonic leaves

**SAM-rooted trajectory analysis and Figure 10C**

This workflow estimates pseudotime across five developing-shoot domains:

`SAM → P1_P2 → P3 → P4 → P5`

The analysis uses the 14 biological replicates while retaining only curated,
non-overlapping spots assigned to these five domains. The biological starting
state is specified as the shoot apical meristem (`SAM`). Monocle 3 then selects
the principal graph node containing the largest number of SAM spots as the root.

The executable script is:

[`01_monocle3_pseudotime_14_samples_Figure_10C.R`](../scripts/R/09_monocle3_pseudotime/01_monocle3_pseudotime_14_samples_Figure_10C.R)

## Inputs

### Individual Seurat objects

The following 14 Seurat v5 objects are required in `data/processed/`:

```text
UL01_seurat_v5.rds    UL02_seurat_v5.rds    UL04_seurat_v5.rds
VR01_seurat_v5.rds    VR02_seurat_v5.rds    VR03_seurat_v5.rds
VR04_seurat_v5.rds    DQ01_seurat_v5.rds    DQ02_seurat_v5.rds
DQ03_seurat_v5.rds    DQ04_seurat_v5.rds    DQ06_seurat_v5.rds
DQ07_seurat_v5.rds    DQ08_seurat_v5.rds
```

Each object must contain an `RNA` assay with one raw `counts` layer.

### Combined Seurat object

```text
data/processed/maize_shoot_14samples_SCT_harmony_seurat_v5.rds
```

For manuscript reproduction, the script transfers the `umap_harmony`
coordinates from this object into Monocle 3. This makes the geometry of Figure
10C correspond to the Harmony UMAP used in the associated Seurat and RNA
velocity panels.

### Curated spot metadata

```text
data/metadata/metadata.csv
```

Required columns are:

- `Barcode`: canonical spot identifier;
- `sample_id`: one of the 14 biological replicates;
- `domains`: anatomical domain; and
- `section_id`: unique physical-section identifier used for alignment.

The script retains only `SAM`, `P1_P2`, `P3`, `P4`, and `P5`. Spots spanning
multiple anatomical domains and spots from the protective coleoptile tissues
are excluded from this developmental trajectory.

## Barcode standardization

The 14 individual objects may use sample-prefixed barcodes. Before combining
the Monocle datasets, the script converts every name to the canonical format:

```text
<10x-barcode>-1_1_<sample-number>
```

The fixed sample-number order is:

```text
UL01=1, UL02=2, UL04=3,
VR01=4, VR02=5, VR03=6, VR04=7,
DQ01=8, DQ02=9, DQ03=10, DQ04=11,
DQ06=12, DQ07=13, DQ08=14
```

The script stops if canonicalization creates duplicates or if no curated spots
match an individual object.

## Analysis workflow

### 1. Construct and combine Monocle datasets

Raw RNA counts, curated spot metadata, and gene identifiers are used to create
one `cell_data_set` per biological replicate. The 14 objects are combined with:

```r
cds <- combine_cds(
  cds_list,
  keep_all_genes = TRUE,
  cell_names_unique = TRUE,
  sample_col_name = "sample_id_from_combine",
  keep_reduced_dims = FALSE
)
```

Raw total UMI counts are calculated for every retained spot. Spots with fewer
than 100 UMIs are removed, reproducing the original `umi_cutoff = 100`
criterion.

### 2. Preprocess and align by physical section

The combined dataset is log-normalized and processed using 100 principal
components:

```r
cds <- preprocess_cds(
  cds,
  method = "PCA",
  num_dim = 100,
  norm_method = "log"
)
```

Batch alignment uses `section_id`, which uniquely identifies each physical
section:

```r
cds <- align_cds(
  cds,
  preprocess_method = "PCA",
  alignment_group = "section_id"
)
```

Repeated positional labels such as `Section1` are not used as the alignment
group because the same label can occur in different capture areas.

### 3. Transfer the Harmony UMAP for Figure 10C

The default setting is:

```r
use_seurat_harmony_umap <- TRUE
seurat_umap_reduction <- "umap_harmony"
```

The script reads the two-dimensional Harmony UMAP from the combined Seurat
object, matches it to retained Monocle barcodes, and installs it as the Monocle
`UMAP` reduction before clustering and graph learning. This reproduces the
coordinate system of the reference Figure 10C.

For a new dataset, an independent Monocle UMAP can instead be calculated from
the aligned PCA coordinates by changing:

```r
use_seurat_harmony_umap <- FALSE
```

In that mode, the script uses cosine distance, 15 neighbors, and
`umap.min_dist = 0.1`.

### 4. Learn a developmental trajectory

Monocle clustering is calculated on the selected UMAP. A single trajectory
graph is then learned across the five developmental domains:

```r
cds <- cluster_cells(cds, reduction_method = "UMAP")

cds <- learn_graph(
  cds,
  use_partition = FALSE,
  close_loop = FALSE
)
```

`use_partition = FALSE` allows one graph to connect the SAM and developing-leaf
regions. `close_loop = FALSE` retains a tree-like trajectory.

### 5. Define the SAM root and calculate pseudotime

The biological root state is manually defined as `SAM`. The script identifies
the principal graph node nearest to the largest number of SAM spots and passes
that node explicitly to `order_cells()`:

```r
cds <- order_cells(
  cds,
  reduction_method = "UMAP",
  root_pr_nodes = sam_root_node
)
```

This root choice is an explicit biological assumption. Reversing or changing
the root changes the interpretation of pseudotime. Pseudotime is a relative
graph distance and is not chronological time.

## Run the workflow

Restart R to obtain a clean session, open the repository root, and run:

```r
setwd("C:/Users/user/Dropbox/Git/maize_shoot_data_process_v2")

source(
  "scripts/R/09_monocle3_pseudotime/01_monocle3_pseudotime_14_samples_Figure_10C.R"
)
```

Do not source `_verify_monocle_runtime.R`. That name referred to a temporary
diagnostic file and is not part of the repository workflow.

## Outputs

### Processed Monocle object

```text
data/processed/monocle3/
└── maize_shoot_SAM_P1_P2_P3_P4_P5_monocle3_pseudotime.rds
```

### Tables

```text
results/tables/09_monocle3_pseudotime/
├── Seurat_to_Monocle_barcode_matching.csv
├── retained_spots_by_sample_and_domain.csv
├── Monocle3_spot_pseudotime_and_metadata.csv
└── Monocle3_pseudotime_parameters.csv
```

`Monocle3_spot_pseudotime_and_metadata.csv` contains the retained barcode,
sample and section metadata, structural domain, total UMI count, pseudotime,
and UMAP coordinates.

### Figures

```text
results/figures/09_monocle3_pseudotime/
├── Figure_10C_monocle3_pseudotime.png
├── Figure_10C_monocle3_pseudotime.pdf
└── Monocle3_trajectory_by_domain.png
```

## Figure 10C

The following image appears after the workflow has completed:

![Monocle 3 pseudotime trajectory rooted in SAM](../results/figures/09_monocle3_pseudotime/Figure_10C_monocle3_pseudotime.png)

Spots are colored from low to high pseudotime. The dark line represents the
principal graph. Principal points are labeled to document the graph structure,
and anatomical-domain labels are positioned at the median UMAP coordinates for
SAM, P1_P2, P3, P4, and P5.

### Anatomical-domain diagnostic

![Monocle 3 trajectory colored by structural domain](../results/figures/09_monocle3_pseudotime/Monocle3_trajectory_by_domain.png)

This diagnostic should be inspected to confirm that the SAM root is located in
the expected anatomical region and that increasing pseudotime broadly follows
the proposed progression toward later leaf primordia.

## Suggested figure legend

**Figure 10C. Pseudotime trajectory of maize embryonic leaf development.**
Monocle 3 trajectory inferred from spatial transcriptomic spots assigned to the
shoot apical meristem (SAM) and successive leaf primordia (P1–P2, P3, P4, and
P5). The Harmony UMAP coordinates from the processed Seurat object were
transferred to Monocle 3 for graph learning and visualization. The trajectory
was rooted in the principal graph node containing the largest number of SAM
spots. Colors indicate increasing pseudotime, the dark line represents the
learned principal graph, and numbered labels identify principal graph points.

## Interpretation and reproducibility notes

- The SAM root is biologically specified rather than inferred without prior
  information.
- Spots from the same section or capture area are subsamples, not independent
  biological replicates.
- The trajectory describes transcriptional continuity and does not establish
  lineage relationships or elapsed developmental time by itself.
- Harmony UMAP transfer is used to reproduce the manuscript coordinate system.
  New datasets should be reprocessed and evaluated independently.
- Principal-point numbers are graph identifiers and are not anatomical-stage
  labels.
- The random seed is recorded as `2026`; numerical results may still vary with
  changes in Monocle 3, UMAP, or graph-learning dependencies.

## Troubleshooting

### S4 object cannot be coerced to a character vector

The script explicitly uses `SeuratObject::Assays()`,
`SeuratObject::Layers()`, and `SeuratObject::LayerData()` to avoid a namespace
collision with Bioconductor assay accessors. Restart R and source the current
repository script if this error appears.

### The combined object lacks `umap_harmony`

Run the step 04 integration workflow first or set:

```r
use_seurat_harmony_umap <- FALSE
```

The latter calculates a new Monocle UMAP and will not reproduce the exact
Figure 10C coordinate system.

### Infinite pseudotime values

Infinite values indicate spots that are not reachable from the selected SAM
root. Confirm that `use_partition = FALSE`, inspect the domain-colored graph,
and verify that the SAM and developing-leaf regions are connected appropriately.

## Session information

The script writes the full R session record to:

[`09_monocle3_pseudotime_sessionInfo.txt`](../results/sessionInfo/09_monocle3_pseudotime_sessionInfo.txt)

```r
sessionInfo()
```
