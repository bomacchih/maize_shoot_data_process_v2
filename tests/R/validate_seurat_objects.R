#!/usr/bin/env Rscript

# Validate the structure of the Seurat objects declared present by the project.
# This test is intentionally separate from the standard-library Python checks
# because reading RDS and Seurat v5 layers requires the configured R environment.

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
})

find_project_root <- function(path = getwd()) {
  path <- normalizePath(path, winslash = "/", mustWork = TRUE)
  repeat {
    if (file.exists(file.path(path, "config", "spaceranger_samples.csv")) &&
        dir.exists(file.path(path, "scripts"))) return(path)
    parent <- dirname(path)
    if (identical(parent, path)) stop("Could not find the project root.")
    path <- parent
  }
}

failures <- character()
passes <- character()

check <- function(condition, message) {
  if (!isTRUE(condition)) failures <<- c(failures, message)
}

project_root <- find_project_root()
manifest <- read.csv(
  file.path(project_root, "config", "spaceranger_samples.csv"),
  stringsAsFactors = FALSE,
  check.names = FALSE
)

expected_samples <- c(
  "UL01", "UL02", "UL04",
  "VR01", "VR02", "VR03", "VR04",
  "DQ01", "DQ02", "DQ03", "DQ04", "DQ06", "DQ07", "DQ08"
)
check(identical(manifest$sample_id, expected_samples),
      "Sample manifest order differs from the fixed 14-sample order.")
check(identical(as.integer(manifest$barcode_suffix), seq_along(expected_samples)),
      "Barcode suffixes are not the fixed 1-14 sequence.")

for (sample_id in expected_samples) {
  path <- file.path(
    project_root, "data", "processed", paste0(sample_id, "_seurat_v5.rds")
  )
  if (!file.exists(path)) {
    failures <- c(failures, paste("Missing individual object:", path))
    next
  }
  object <- tryCatch(readRDS(path), error = function(error) error)
  if (inherits(object, "error")) {
    failures <- c(failures, paste(sample_id, "could not be read:", object$message))
    next
  }
  check(inherits(object, "Seurat"), paste(sample_id, "is not a Seurat object."))
  assays <- Assays(object)
  check(any(c("RNA", "Spatial") %in% assays),
        paste(sample_id, "has neither an RNA nor Spatial assay."))
  check(ncol(object) > 0L, paste(sample_id, "contains no spots."))
  check(nrow(object) > 0L, paste(sample_id, "contains no genes."))
  passes <- c(passes, paste(sample_id, "read successfully"))
}

integrated_path <- file.path(
  project_root, "data", "processed",
  "maize_shoot_14samples_SCT_harmony_seurat_v5.rds"
)
if (!file.exists(integrated_path)) {
  failures <- c(failures, paste("Missing integrated object:", integrated_path))
} else {
  integrated <- tryCatch(readRDS(integrated_path), error = function(error) error)
  if (inherits(integrated, "error")) {
    failures <- c(failures, paste("Integrated object could not be read:", integrated$message))
  } else {
    check(inherits(integrated, "Seurat"), "Integrated object is not a Seurat object.")
    check(all(c("RNA", "SCT") %in% Assays(integrated)),
          "Integrated object does not contain both RNA and SCT assays.")
    check(all(c("pca", "harmony", "umap_pca", "umap_harmony") %in%
                Reductions(integrated)),
          "Integrated object is missing an expected PCA/Harmony/UMAP reduction.")
    required_metadata <- c(
      "Barcode", "sample_id", "section_id", "domains", "harmony_clusters",
      "percent.mito", "percent.pltd"
    )
    check(all(required_metadata %in% colnames(integrated[[]])),
          "Integrated object is missing required metadata fields.")
    if ("Barcode" %in% colnames(integrated[[]])) {
      check(identical(as.character(integrated$Barcode), colnames(integrated)),
            "Integrated Barcode metadata is not aligned to object column names.")
    }
    if ("sample_id" %in% colnames(integrated[[]])) {
      check(setequal(unique(as.character(integrated$sample_id)), expected_samples),
            "Integrated object does not contain exactly the 14 configured samples.")
    }
    passes <- c(passes, "integrated object read successfully")
  }
}

for (message in passes) cat("[PASS]", message, "\n")
if (length(failures)) {
  for (message in failures) cat("[FAIL]", message, "\n")
  cat("\n", length(failures), "Seurat validation failure(s).\n")
  quit(status = 1L)
}

cat("\nAll", length(expected_samples),
    "individual objects and the integrated Seurat object passed structural validation.\n")

