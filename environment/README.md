# Reproducible software environments

The repository shares environment specifications rather than installed library
or environment directories. The latter contain machine-specific paths and
binaries and must not be copied between computers.

## R environment

The locked environment uses R 4.6.1. From the repository root, install `renv`
and restore the project library:

```r
install.packages("renv")
renv::restore()
renv::status()
```

`DESCRIPTION` declares the direct R dependencies, while `renv.lock` records
their exact versions, sources, repositories, and transitive dependencies. The
project-local `renv/library/` directory is intentionally ignored by Git.

After deliberately changing and testing an R package, update the lockfile with:

```r
renv::snapshot()
```

Audit the active R environment with:

```bash
Rscript tests/R/validate_R_environment.R
```

## Python environment

The tested Python/scVelo environment uses Python 3.11 and Conda Forge. Create
the portable environment with:

```bash
conda env create --file environment/python/environment.yml
conda activate maize-shoot-v2
python scripts/python/08_RNA_velocity/01_scvelo_dynamical_RNA_velocity.py --help
```

For an exact reconstruction on Windows x86-64, use the fully resolved platform
environment instead:

```bash
conda env create --file environment/python/environment-win-64.yml
```

`environment.yml` is the cross-platform sharing specification.
`environment-win-64.yml` records exact package versions and Conda build strings
for the validated Windows environment. Do not commit the Conda environment
directory itself.

The Python workflow records the versions actually used during an analysis in
`results/tables/08_RNA_velocity/python_package_versions.json`.

## SCINA

SCINA 1.2.0 was removed from the active CRAN repository on 2026-05-08, so
`install.packages("SCINA")` no longer works against the normal CRAN package
index. The source package remains available from the official CRAN archive.

Install the pinned archived release with:

```bash
Rscript environment/install_SCINA.R
```

The installer uses the official HTTPS archive URL, installs into the first
active R library, verifies version 1.2.0, and does nothing when the correct
version is already installed.
