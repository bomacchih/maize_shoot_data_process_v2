# R environment notes

Validated on 2026-08-30 with R 4.6.1.

Run the environment audit from the repository root:

```bash
Rscript tests/R/validate_R_environment.R
```

The current inventory contains 26 required package-version checks. The earlier
audit matched all 26 required packages with zero version
warnings, and zero missing required packages. This includes Shiny 1.14.0 for
the interactive Monocle 3 root-principal-node selection workflow,
`sctransform` 0.4.3 for SCT model handling, and `presto` 1.0.0 for sparse AUC
marker scoring. Rerun the command above after changing the active R library;
the report is environment-specific rather than a permanent installation claim.

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
