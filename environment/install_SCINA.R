#!/usr/bin/env Rscript

scina_version <- "1.2.0"
scina_url <- paste0(
  "https://cran.r-project.org/src/contrib/Archive/SCINA/",
  "SCINA_", scina_version, ".tar.gz"
)
target_library <- .libPaths()[1L]

dependencies <- c("MASS", "gplots")
missing_dependencies <- dependencies[
  !vapply(dependencies, requireNamespace, logical(1L), quietly = TRUE)
]
if (length(missing_dependencies)) {
  stop(
    "Install the required package(s) first: ",
    paste(missing_dependencies, collapse = ", ")
  )
}

installed_correct_version <- requireNamespace("SCINA", quietly = TRUE) &&
  identical(as.character(packageVersion("SCINA")), scina_version)

if (!installed_correct_version) {
  dir.create(target_library, recursive = TRUE, showWarnings = FALSE)
  install.packages(
    scina_url,
    repos = NULL,
    type = "source",
    lib = target_library
  )
}

if (!requireNamespace("SCINA", quietly = TRUE)) {
  stop("SCINA could not be loaded after installation.")
}
if (!identical(as.character(packageVersion("SCINA")), scina_version)) {
  stop(
    "Expected SCINA ", scina_version, " but found ",
    as.character(packageVersion("SCINA")), "."
  )
}

message(
  "SCINA ", scina_version, " is available at ",
  find.package("SCINA")
)
