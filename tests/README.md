# Validation tests

Run commands from the repository root.

## Configuration and input validation

The main validator uses only the Python standard library:

```bash
python scripts/python/validate_project_config.py
```

This checks configuration syntax and ranges, the fixed sample order, barcode
suffixes, metadata integrity, files declared present, and consistency between
the configuration snapshot and important script assignments.

Warnings identify unresolved scientific settings or unavailable end-to-end
inputs. To make warnings fail a release check:

```bash
python scripts/python/validate_project_config.py --strict
```

Write a machine-readable report with:

```bash
python scripts/python/validate_project_config.py \
  --json-report results/logs/config_validation_report.json
```

## Python unit tests

```bash
python -m unittest discover -s tests/python -p "test_*.py"
```

Check the active Python analysis environment separately:

```bash
python tests/python/validate_python_environment.py
```

This reports missing packages and also warns that the current Python
requirements are not yet pinned to exact versions.

## Cell-type marker validation

The committed SCINA marker CSV is validated automatically by the main Python
validator. To also verify the original curated RDS against the expression
features in both the combined scRNA reference and a representative Visium
object, run:

```bash
Rscript tests/R/validate_marker_list_for_mapping.R \
  <path-to-marker_list2.rds> \
  data/processed/sc_merged_filter_SCT2_inte.rds \
  data/processed/UL01_seurat_v5.rds
```

The check fails for blank or duplicated markers, unstable maize gene IDs,
markers assigned to more than one cell type, or markers absent from either
feature space.

Compare the committed marker file directly with Supplementary Table 9 using
only the Python standard library:

```bash
python scripts/python/compare_supplementary_table9_markers.py \
  <path-to-Supplementary_Tables_20251104.xlsx> \
  data/metadata/scRNA_reference/SCINA_marker_table.csv
```

This confirms that every retained marker is an exclusive Supplementary Table
9 assignment and that its source order is preserved after filtering.

## Seurat object validation

Run this only after activating the documented R environment:

```bash
Rscript tests/R/validate_seurat_objects.R
```

It loads the 14 individual objects and the integrated object, then checks their
assays, reductions, required metadata, sample coverage, and barcode alignment.

## R environment validation

After activating the intended R environment, compare installed packages with
the audited inventory:

```bash
Rscript tests/R/validate_R_environment.R
```

Missing required packages fail the check. Version differences are reported as
warnings until the final `renv.lock` environment is created.

## Interpretation

- `PASS` means an objective check succeeded.
- `WARN` means the repository can be inspected but a scientific choice or
  external input remains unresolved.
- `FAIL` means a declared invariant is broken and should be fixed before an
  analysis or release.
