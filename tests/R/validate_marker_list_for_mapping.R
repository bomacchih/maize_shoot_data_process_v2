#!/usr/bin/env Rscript

# Validate a marker-list RDS against the scRNA and Visium feature spaces.

arguments <- commandArgs(trailingOnly = TRUE)
if (length(arguments) != 3L) {
  stop(
    "Usage: Rscript validate_marker_list_for_mapping.R ",
    "<marker_list.rds> <scRNA_reference.rds> <Visium_object.rds>"
  )
}

marker_list <- readRDS(arguments[[1L]])
scrna_reference <- readRDS(arguments[[2L]])
visium_object <- readRDS(arguments[[3L]])

if (!is.list(marker_list) || is.null(names(marker_list))) {
  stop("The marker RDS must contain a named list.")
}
marker_list <- lapply(marker_list, function(markers) {
  unique(trimws(as.character(markers)))
})
markers <- unique(unlist(marker_list, use.names = FALSE))
scrna_genes <- rownames(scrna_reference)
visium_genes <- rownames(visium_object)

if (!length(scrna_genes)) stop("The scRNA reference has no expression features.")
if (!length(visium_genes)) stop("The Visium object has no expression features.")

missing_scrna <- setdiff(markers, scrna_genes)
missing_visium <- setdiff(markers, visium_genes)

cat("cell_types=", length(marker_list), "\n", sep = "")
cat("unique_markers=", length(markers), "\n", sep = "")
cat("scRNA_overlap=", length(markers) - length(missing_scrna), "/",
    length(markers), "\n", sep = "")
cat("Visium_overlap=", length(markers) - length(missing_visium), "/",
    length(markers), "\n", sep = "")

if (length(missing_scrna)) {
  stop(length(missing_scrna), " marker(s) are absent from the scRNA reference.")
}
if (length(missing_visium)) {
  stop(length(missing_visium), " marker(s) are absent from the Visium object.")
}

cat("Marker-list feature-space validation passed.\n")

