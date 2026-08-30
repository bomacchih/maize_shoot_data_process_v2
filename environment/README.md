# R environment notes

Validated on 2026-08-30 with R 4.6.1.

Run the environment audit from the repository root:

```bash
Rscript tests/R/validate_R_environment.R
```

The current audit result is 24 exact package-version matches, zero version
warnings, and zero missing required packages. This includes Shiny 1.14.0 for
the interactive Monocle 3 root-principal-node selection workflow.

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
