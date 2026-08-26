# Maize shoot spatial transcriptomics data-processing pipeline

This repository contains reproducible R, Python, and shell workflows for processing maize developing-shoot spatial transcriptomics data and integrating it with a maize scRNA-seq reference. The workflows support the revised Bio-protocol manuscript associated with the original Plant Biotechnology Journal study.

## Main workflow

1. Build the maize Space Ranger reference and process Visium libraries.
2. Load individual Space Ranger outputs into STUtility or Seurat v5.
3. Add curated spot metadata.
4. Perform quality control, SCTransform normalization, and Harmony integration across 14 capture areas.
5. Identify markers and assign anatomical tissue supergroups.
6. Generate structural-domain pseudobulk profiles.
7. Cluster developmental expression trends and perform GO enrichment.
8. Estimate RNA velocity with Velocyto and scVelo.
9. Estimate pseudotime with Monocle 3.
10. Prepare and annotate the maize scRNA-seq reference and map cell-type information to Visium spots with Seurat and SPOTlight.

## Repository structure

```text
config/       Configuration files
data/         Dataset and software inventories, metadata, and small references
docs/         Workflow documentation
environment/  Environment specifications
results/      Generated figures, tables, and session information
scripts/      R, Python, and shell pipelines
tests/        Validation scripts
workflow/     Workflow-level orchestration files
```

## Data availability

Large raw and processed datasets are intentionally excluded from Git. See:

- [`data/DATASETS.md`](data/DATASETS.md) for required datasets, expected locations, current availability, and missing inputs.
- [`data/dataset_inventory.csv`](data/dataset_inventory.csv) for a machine-readable inventory.

## Software requirements

- [`data/SOFTWARE_REQUIREMENTS.md`](data/SOFTWARE_REQUIREMENTS.md) summarizes the R, Python, and external command-line requirements.
- [`data/R_packages.csv`](data/R_packages.csv) records the R packages and versions available in the audit environment.
- [`data/requirements-python.txt`](data/requirements-python.txt) lists Python dependencies for the RNA-velocity workflow.

## Starting the analysis

Run scripts from the repository root so that project-relative paths resolve correctly. Begin with the dataset inventory to determine whether the raw or processed entry point is available. Each completed workflow writes figures and tables under `results/` and records runtime information where implemented.

## Related publication

The protocol was used in the following study:

Wu et al. *Plant Biotechnology Journal* (2026). DOI: [10.1111/pbi.70515](https://doi.org/10.1111/pbi.70515).

## Repository scope

The code is under active manuscript revision. Large Seurat RDS objects, Space Ranger outputs, FASTQ files, loom files, and other sequencing-derived binaries must be obtained from the associated data archive and placed according to `data/DATASETS.md`.
