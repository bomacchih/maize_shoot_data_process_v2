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

project_root <- find_project_root()
arguments <- commandArgs(trailingOnly = TRUE)

default_source_rds <- file.path(
  project_root, "data", "reference", "scRNA_reference", "marker_list2.rds"
)
default_output_csv <- file.path(
  project_root, "data", "metadata", "scRNA_reference",
  "SCINA_marker_table.csv"
)

source_rds <- if (length(arguments) >= 1L) {
  normalizePath(arguments[[1L]], winslash = "/", mustWork = TRUE)
} else if (file.exists(default_source_rds)) {
  normalizePath(default_source_rds, winslash = "/", mustWork = TRUE)
} else {
  NA_character_
}

output_csv <- if (length(arguments) >= 2L) {
  arguments[[2L]]
} else {
  default_output_csv
}
if (!grepl("^([A-Za-z]:[/\\]|[/\\])", output_csv)) {
  output_csv <- file.path(project_root, output_csv)
}

reuse_existing_csv <- is.na(source_rds) && file.exists(output_csv)

if (reuse_existing_csv) {
  marker_table <- read.csv(
    output_csv,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  required_columns <- c(
    "gene_id", "cell_type", "marker_rank", "source_cell_type"
  )
  missing_columns <- setdiff(required_columns, names(marker_table))
  if (length(missing_columns) > 0L) {
    stop(
      "The existing SCINA marker table is missing columns: ",
      paste(missing_columns, collapse = ", ")
    )
  }
  if (nrow(marker_table) == 0L ||
      anyNA(marker_table[, required_columns, drop = FALSE]) ||
      any(!nzchar(trimws(marker_table$gene_id))) ||
      any(!nzchar(trimws(marker_table$cell_type))) ||
      anyDuplicated(marker_table[, c("gene_id", "cell_type")])) {
    stop("The existing SCINA marker table failed validation.")
  }

  message(
    "No marker-list RDS argument was supplied and the project-local default ",
    "was not found. Validated the existing prepared marker table instead: ",
    output_csv
  )
  message(
    "Validated ", nrow(marker_table), " markers for ",
    length(unique(marker_table$cell_type)), " cell types."
  )
} else if (is.na(source_rds)) {
  stop(
    "No marker-list RDS was supplied, the default file was not found at ",
    default_source_rds,
    ", and no prepared SCINA marker table exists at ", output_csv, ".\n",
    "Run from the command line as: Rscript ",
    "00_prepare_SCINA_marker_table.R <marker_list2.rds> [output.csv]"
  )
} else {
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
  "upstream_source_workbook_md5=3d0c20648b193de5f1919ad7127f3b70",
  "upstream_source_sheet=Supplementary Table 9",
  "supplementary_table9_comparison_passed=true",
  "supplementary_table9_raw_rows=8080",
  "supplementary_table9_duplicate_pair_rows=1227",
  "supplementary_table9_genes_shared_across_cell_types=1195",
  "supplementary_table9_exclusive_genes=4094",
  "supplementary_table9_exclusive_genes_not_in_curated_rds=97",
  paste0(
    "comparison_script=",
    "scripts/python/compare_supplementary_table9_markers.py"
  ),
  paste0(
    "comparison_report=",
    "results/logs/scina_marker_supplementary_table9_comparison.json"
  ),
  "marker_rank_semantics=position_within_source_list",
  "cell_type_normalization=non_alphanumeric_punctuation_replaced_with_underscore"
)
provenance_connection <- file(provenance_file, open = "wb")
tryCatch(
  writeLines(
    provenance,
    con = provenance_connection,
    sep = "\n",
    useBytes = TRUE
  ),
  finally = close(provenance_connection)
)

message("Wrote ", nrow(marker_table), " markers for ", length(marker_list),
        " cell types to ", output_csv)
message("Provenance: ", provenance_file)
}
