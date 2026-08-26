# R and Python package requirements

Audit date: 2026-08-26

This list was extracted from the scripts in this repository. The installed R versions are the versions available in the environment used for the code audit; they are not strict version pins unless explicitly stated. Every completed analysis should retain its own `sessionInfo()` or Python version record.

## R environment

- Audit environment: **R 4.6.1**
- Pipeline target: R 4.x with Seurat v5-compatible packages
- Machine-readable list: `data/R_packages.csv`

### Required R packages

| Package | Audit version | Main use |
|---|---:|---|
| Seurat | 5.5.1 | Visium and scRNA-seq objects, SCTransform, integration, clustering, marker analysis, and plotting |
| SeuratObject | 5.4.0 | Seurat v5 assays, layers, identities, reductions, and metadata access |
| STutility | 1.1.1 | Alternative loading and handling of Space Ranger spatial data |
| Matrix | 1.7-6 | Sparse-matrix operations |
| ggplot2 | 4.0.3 | Figures |
| patchwork | 1.3.2 | Multi-panel figures |
| harmony | 2.0.5 | Batch integration through Harmony |
| dplyr | 1.2.1 | Table and metadata manipulation |
| tidyr | 1.3.2 | Long/wide reshaping for developmental-trend analysis |
| edgeR | 4.10.1 | TMM normalization and log2 CPM pseudobulk profiles |
| monocle3 | 1.4.27 | Pseudotime and trajectory analysis |
| SingleCellExperiment | 1.34.0 | scRNA/Monocle/SPOTlight data exchange |
| SummarizedExperiment | 1.42.0 | Assay and column-metadata access |
| S4Vectors | 0.50.1 | Bioconductor metadata containers used during mapping |
| igraph | 2.3.3 | Monocle principal-graph node handling |
| scDblFinder | 1.26.7 | Doublet detection in scRNA-seq reference libraries |
| glmGamPoi | 1.24.0 | SCTransform v2 fitting method |
| SpatialExperiment | 1.22.0 | Spatial data containers for SPOTlight |
| scuttle | 1.22.0 | scRNA-seq log normalization for SPOTlight reference preparation |
| scran | 1.40.0 | HVG and marker scoring for SPOTlight |
| SPOTlight | 1.16.0 | Visium spot cell-type deconvolution and scatter-pie plots |
| SCINA | **1.2.0** | Cell-type annotation; the script explicitly validates against version 1.2.0 |
| scales | 1.4.0 | Continuous color-scale bounds in SPOTlight UMAP plots |

### Optional R packages

| Package | Audit version | Behavior when absent |
|---|---:|---|
| ggrepel | 0.9.8 | Labels use simpler placement when unavailable |
| GO.db | Not installed | GO descriptions must come from the supplied reference CSV/mapping |
| AnnotationDbi | 1.74.0 | Used with `GO.db` to retrieve GO term descriptions |

Base and recommended R packages called explicitly through namespaces—such as `stats`, `utils`, `grDevices`, and `grid`—ship with R and do not require separate installation.

## Python environment

- Recommended interpreter: **Python 3.10–3.12**
- Pip-compatible list: `data/requirements-python.txt`
- The RNA-velocity script writes detected versions to `results/tables/08_RNA_velocity/python_package_versions.json` after a successful run.

### Required Python packages

| Package | Import name | Main use |
|---|---|---|
| anndata | `anndata` | AnnData containers and H5AD output |
| matplotlib | `matplotlib` | RNA-velocity figures |
| numpy | `numpy` | Numerical arrays |
| pandas | `pandas` | Metadata and exported tables |
| scanpy | `scanpy` | Preprocessing, PCA, neighbors, and UMAP |
| scvelo | `scvelo` | Dynamical RNA-velocity modeling |
| scipy | `scipy` | Sparse matrices and numerical routines |
| loompy | `loompy` | Reading Velocyto loom files |
| h5py | `h5py` | HDF5/H5AD input-output support |
| umap-learn | `umap` | UMAP calculation used through Scanpy |
| scikit-learn | `sklearn` | Neighbor and dimensionality-reduction support used through Scanpy |

### Optional Python package

| Package | When required |
|---|---|
| harmonypy | Only when running the RNA-velocity script with `--use-harmony` |

## External command-line software

| Software | Role |
|---|---|
| 10x Genomics Space Ranger | Reference construction and Visium `count` processing |
| Velocyto.py 0.17 series | Generation of spliced/unspliced loom files before scVelo analysis |

These command-line tools are not installed through the R or Python package lists. Their executable versions and reference-build parameters should be recorded separately.

## Reproducibility recommendation

Do not treat the audit versions as permanent requirements. After the full pipeline runs successfully, record:

1. R and all R package versions with `sessionInfo()`.
2. Python and package versions using the JSON file produced by the RNA-velocity script.
3. Space Ranger and Velocyto versions in the shell-processing log.
4. A lock file or environment specification—such as `renv.lock` for R and a pinned Python requirements/environment file—for the final release.
