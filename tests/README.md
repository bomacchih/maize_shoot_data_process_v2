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
