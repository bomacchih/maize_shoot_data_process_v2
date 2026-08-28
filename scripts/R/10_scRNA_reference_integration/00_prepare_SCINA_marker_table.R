#!/usr/bin/env Rscript

# Convert the curated marker_list2.rds named list into the CSV consumed by the
# SCINA annotation workflow. The order of genes within each source list is
# preserved as marker_rank; it is not claimed to be an effect-size ranking.

find_project_root <- function(start = getwd()) {
  current <- normalizePath(start, winslash = "/", mustWork = TRUE)
  repeat {
    if (dir.exists(file.path(current, "config")) &&
        dir.exists(file.path(current, "scripts"))) {
      return(current)
    }
    parent <- dirname(current)
    if (identical(parent, current)) {
      stop("Could not locate the repository root from: ", start)
    }
    current <- parent
  }
}

normalize_cell_type <- function(value) {
  normalized <- gsub("[^A-Za-z0-9_]+", "_", value)
  normalized <- gsub("_+", "_", normalized)
  gsub("^_|_$", "", normalized)
}

arguments <- commandArgs(trailingOnly = TRUE)
if (!length(arguments)) {
  stop(
    "Usage: Rscript 00_prepare_SCINA_marker_table.R ",
    "<marker_list2.rds> [output.csv]"
  )
}

project_root <- find_project_root()
source_rds <- normalizePath(arguments[[1L]], winslash = "/", mustWork = TRUE)
output_csv <- if (length(arguments) >= 2L) {
  arguments[[2L]]
} else {
  file.path(
    project_root, "data", "metadata", "scRNA_reference",
    "SCINA_marker_table.csv"
  )
}
if (!grepl("^([A-Za-z]:[/\\]|[/\\])", output_csv)) {
  output_csv <- file.path(project_root, output_csv)
}

marker_list <- readRDS(source_rds)
if (!is.list(marker_list) || !length(marker_list)) {
  stop("The marker RDS must contain a non-empty list.")
}
if (is.null(names(marker_list)) || anyNA(names(marker_list)) ||
    any(!nzchar(trimws(names(marker_list))))) {
  stop("Every marker-list element must have a non-empty cell-type name.")
}
if (anyDuplicated(names(marker_list))) {
  stop("Marker-list cell-type names must be unique.")
}

marker_list <- lapply(marker_list, function(markers) {
  markers <- trimws(as.character(markers))
  if (anyNA(markers) || any(!nzchar(markers))) {
    stop("Marker lists cannot contain missing or blank gene identifiers.")
  }
  if (anyDuplicated(markers)) {
    stop("Marker lists cannot contain duplicate gene identifiers.")
  }
  if (any(!grepl("^Zm[0-9]+[A-Za-z]+[0-9]+$", markers))) {
    stop("Every marker must use a stable maize gene identifier beginning with Zm.")
  }
  markers
})

normalized_cell_types <- vapply(
  names(marker_list), normalize_cell_type, character(1L)
)
if (anyDuplicated(normalized_cell_types)) {
  stop("Cell-type normalization produced duplicate labels.")
}

gene_membership <- table(unlist(marker_list, use.names = FALSE))
shared_gene_count <- sum(gene_membership > 1L)
if (shared_gene_count > 0L) {
  stop(
    shared_gene_count,
    " gene(s) occur in more than one cell-type list; resolve overlap before conversion."
  )
}

marker_table <- do.call(
  rbind,
  lapply(seq_along(marker_list), function(index) {
    genes <- marker_list[[index]]
    data.frame(
      gene_id = genes,
      cell_type = normalized_cell_types[[index]],
      marker_rank = seq_along(genes),
      source_cell_type = names(marker_list)[[index]],
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  })
)
rownames(marker_table) <- NULL

dir.create(dirname(output_csv), recursive = TRUE, showWarnings = FALSE)
write.csv(marker_table, output_csv, row.names = FALSE, quote = TRUE)

provenance_file <- file.path(dirname(output_csv), "SCINA_marker_table_provenance.txt")
provenance <- c(
  paste0("source_file=", basename(source_rds)),
  paste0("source_md5=", unname(tools::md5sum(source_rds))),
  paste0("source_cell_types=", length(marker_list)),
  paste0("marker_assignments=", nrow(marker_table)),
  paste0("unique_genes=", length(gene_membership)),
  paste0("shared_genes=", shared_gene_count),
  "upstream_source_workbook=Supplementary_Tables_20251104.xlsx",
  "upstream_source_sheet=Supplementary Table 9",
  paste0(
    "comparison_script=",
    "scripts/python/compare_supplementary_table9_markers.py"
  ),
  "marker_rank_semantics=position_within_source_list",
  "cell_type_normalization=non_alphanumeric_punctuation_replaced_with_underscore"
)
writeLines(provenance, provenance_file, useBytes = TRUE)

message("Wrote ", nrow(marker_table), " markers for ", length(marker_list),
        " cell types to ", output_csv)
message("Provenance: ", provenance_file)
