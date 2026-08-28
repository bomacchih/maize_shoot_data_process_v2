#!/usr/bin/env Rscript

# Check the active R environment against data/R_packages.csv without loading
# every package. Run this after activating the intended analysis environment.

find_project_root <- function(path = getwd()) {
  path <- normalizePath(path, winslash = "/", mustWork = TRUE)
  repeat {
    if (file.exists(file.path(path, "data", "R_packages.csv")) &&
        dir.exists(file.path(path, "scripts"))) return(path)
    parent <- dirname(path)
    if (identical(parent, path)) stop("Could not find the project root.")
    path <- parent
  }
}

project_root <- find_project_root()
inventory <- read.csv(
  file.path(project_root, "data", "R_packages.csv"),
  stringsAsFactors = FALSE,
  check.names = FALSE
)

required_columns <- c("package", "requirement", "audit_version")
missing_columns <- setdiff(required_columns, colnames(inventory))
if (length(missing_columns)) {
  stop("R_packages.csv is missing: ", paste(missing_columns, collapse = ", "))
}

required <- startsWith(inventory$requirement, "required")
failures <- character()
warnings <- character()
passes <- character()

for (index in which(required)) {
  package_name <- inventory$package[index]
  expected_version <- inventory$audit_version[index]
  if (!requireNamespace(package_name, quietly = TRUE)) {
    failures <- c(failures, paste(package_name, "is not installed"))
    next
  }
  observed_version <- as.character(packageVersion(package_name))
  versions_match <- expected_version == "" || expected_version == "NOT_INSTALLED" ||
    isTRUE(package_version(observed_version) == package_version(expected_version))
  if (!versions_match) {
    warnings <- c(
      warnings,
      paste(package_name, "is", observed_version, "but audit version is", expected_version)
    )
  } else {
    passes <- c(passes, paste(package_name, observed_version))
  }
}

cat("[INFO] R", as.character(getRversion()), "\n")
for (message in passes) cat("[PASS]", message, "\n")
for (message in warnings) cat("[WARN]", message, "\n")
for (message in failures) cat("[FAIL]", message, "\n")

cat(
  "\nSummary:", length(passes), "matching,",
  length(warnings), "version warning(s),",
  length(failures), "missing required package(s).\n"
)

if (length(failures)) quit(status = 1L)
