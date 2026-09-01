#!/usr/bin/env Rscript

# Export domain-level RNA expression matrices for the TO-GCN workflow.
#
# Input:
#   data/processed/XGE202122_S5_subset_embleaf_harmony_join.rds
#
# Grouping variable:
#   domains, ordered as SAM -> P1_P2 -> P3 -> P4 -> P5
#
# Outputs:
#   results/tables/11_TO_GCN/input/domain_spot_counts.csv
#   results/tables/11_TO_GCN/input/
#       RNA_AggregateExpression_sum_counts_by_domain.csv.gz
#   results/tables/11_TO_GCN/input/
#       RNA_AverageExpression_mean_counts_by_domain.csv.gz
#   results/tables/11_TO_GCN/input/
#       RNA_Aggregate_vs_Average_validation.csv
#   results/sessionInfo/11_TO_GCN_domain_expression_export_sessionInfo.txt
#
# AggregateExpression() returns summed raw counts for each domain.
# AverageExpression(layer = "counts") returns mean raw counts per spot without
# exponentiating the values. The script confirms that:
#
#   aggregate count / number of spots == average count
#
# Run from the repository root:
#   source("scripts/R/11_TO_GCN/00_export_RNA_expression_by_domain_Seurat_v5.R")

suppressPackageStartupMessages({
    library(Seurat)
    library(SeuratObject)
    library(Matrix)
})

find_project_root <- function(path = getwd()) {
    path <- normalizePath(path, winslash = "/", mustWork = TRUE)
    repeat {
        if (dir.exists(file.path(path, "data")) &&
            dir.exists(file.path(path, "scripts"))) {
            return(path)
        }
        parent <- dirname(path)
        if (identical(parent, path)) {
            stop("Project root not found. Start R in the repository or a subdirectory.")
        }
        path <- parent
    }
}

detect_column <- function(column_names, candidates) {
    normalized <- gsub("[^a-z0-9]+", "_", tolower(trimws(column_names)))
    for (candidate in candidates) {
        hit <- which(normalized == candidate)
        if (length(hit) > 0L) return(column_names[hit[1L]])
    }
    NA_character_
}

get_raw_rna_counts <- function(object) {
    if (!"RNA" %in% Assays(object)) {
        stop("The input object does not contain an RNA assay.")
    }

    layer_names <- Layers(object[["RNA"]])
    count_layers <- grep("^counts($|\\.)", layer_names, value = TRUE)
    if (length(count_layers) == 0L) {
        stop(
            "No raw RNA count layer was found. Available layers: ",
            paste(layer_names, collapse = ", ")
        )
    }

    count_matrices <- lapply(count_layers, function(layer_name) {
        matrix <- LayerData(object, assay = "RNA", layer = layer_name)
        if (!inherits(matrix, "sparseMatrix")) matrix <- as(matrix, "dgCMatrix")
        matrix
    })

    if (length(count_matrices) == 1L) return(count_matrices[[1L]])

    common_features <- Reduce(intersect, lapply(count_matrices, rownames))
    if (length(common_features) == 0L) {
        stop("The RNA count layers have no common features.")
    }
    count_matrices <- lapply(
        count_matrices,
        function(matrix) matrix[common_features, , drop = FALSE]
    )
    combined <- do.call(cbind, count_matrices)
    if (anyDuplicated(colnames(combined))) {
        stop("Duplicated spot names were found while combining RNA count layers.")
    }
    combined
}

write_expression_matrix <- function(expression_matrix, path) {
    output <- data.frame(
        gene_id = rownames(expression_matrix),
        as.data.frame(as.matrix(expression_matrix), check.names = FALSE),
        check.names = FALSE
    )
    connection <- gzfile(path, open = "wt")
    on.exit(close(connection), add = TRUE)
    write.csv(output, connection, row.names = FALSE, quote = FALSE)
}

project_root <- find_project_root()
input_rds <- file.path(
    project_root,
    "data", "processed", "XGE202122_S5_subset_embleaf_harmony_join.rds"
)
output_dir <- file.path(
    project_root, "results", "tables", "11_TO_GCN", "input"
)
session_dir <- file.path(project_root, "results", "sessionInfo")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(session_dir, recursive = TRUE, showWarnings = FALSE)

domain_levels <- c("SAM", "P1_P2", "P3", "P4", "P5")
# Seurat replaces underscores in identity-class names with dashes. Use these
# safe temporary labels during aggregation, then restore the manuscript labels
# in the exported matrices.
seurat_domain_levels <- gsub("_", "-", domain_levels, fixed = TRUE)
domain_column_candidates <- c("domains", "domain", "structural_domain")

if (!file.exists(input_rds)) stop("Input RDS not found: ", input_rds)
embryonic_leaf <- readRDS(input_rds)
stopifnot(inherits(embryonic_leaf, "Seurat"))

raw_counts <- get_raw_rna_counts(embryonic_leaf)
metadata <- embryonic_leaf[[]]
domain_column <- detect_column(
    colnames(metadata), domain_column_candidates
)
if (is.na(domain_column)) {
    stop(
        "Could not identify structural-domain metadata. Available columns: ",
        paste(colnames(metadata), collapse = ", ")
    )
}

missing_metadata_spots <- setdiff(colnames(raw_counts), rownames(metadata))
if (length(missing_metadata_spots) > 0L) {
    stop(
        length(missing_metadata_spots),
        " RNA-count spots are absent from the Seurat metadata."
    )
}

metadata <- metadata[colnames(raw_counts), , drop = FALSE]
domain_values <- trimws(as.character(metadata[[domain_column]]))
keep_spots <- !is.na(domain_values) & domain_values %in% domain_levels
raw_counts <- raw_counts[, keep_spots, drop = FALSE]
metadata <- metadata[keep_spots, , drop = FALSE]
metadata$domains <- factor(domain_values[keep_spots], levels = domain_levels)

if (ncol(raw_counts) == 0L) stop("No SAM–P5 spots remained for export.")
if (!identical(colnames(raw_counts), rownames(metadata))) {
    stop("RNA-count and metadata spot orders are inconsistent.")
}

domain_spot_counts <- table(metadata$domains)
domain_spot_counts <- as.integer(domain_spot_counts[domain_levels])
names(domain_spot_counts) <- domain_levels
if (anyNA(domain_spot_counts) || any(domain_spot_counts == 0L)) {
    stop(
        "At least one required domain has no retained spots: ",
        paste(domain_levels[is.na(domain_spot_counts) | domain_spot_counts == 0L],
              collapse = ", ")
    )
}

# Rebuild a compact RNA-only object so both Seurat functions operate on one
# joined count layer even when the deposited Seurat v5 object stores multiple
# sample-specific RNA count layers.
domain_object <- CreateSeuratObject(
    counts = raw_counts,
    assay = "RNA",
    meta.data = metadata
)
domain_object$domains_for_expression <- factor(
    gsub("_", "-", as.character(domain_object$domains), fixed = TRUE),
    levels = seurat_domain_levels
)

aggregate_expression <- AggregateExpression(
    object = domain_object,
    assays = "RNA",
    group.by = "domains_for_expression",
    return.seurat = FALSE,
    verbose = TRUE
)$RNA

average_expression <- AverageExpression(
    object = domain_object,
    assays = "RNA",
    group.by = "domains_for_expression",
    layer = "counts",
    return.seurat = FALSE,
    verbose = TRUE
)$RNA

missing_aggregate_domains <- setdiff(
    seurat_domain_levels, colnames(aggregate_expression)
)
missing_average_domains <- setdiff(
    seurat_domain_levels, colnames(average_expression)
)
if (length(missing_aggregate_domains) || length(missing_average_domains)) {
    stop(
        "Seurat expression output is missing one or more required domains. ",
        "AggregateExpression returned: ",
        paste(colnames(aggregate_expression), collapse = ", "),
        "; AverageExpression returned: ",
        paste(colnames(average_expression), collapse = ", ")
    )
}
aggregate_expression <- aggregate_expression[
    , seurat_domain_levels, drop = FALSE
]
average_expression <- average_expression[
    , seurat_domain_levels, drop = FALSE
]
colnames(aggregate_expression) <- domain_levels
colnames(average_expression) <- domain_levels

if (!identical(rownames(aggregate_expression), rownames(average_expression))) {
    stop("AggregateExpression and AverageExpression returned different genes.")
}

mean_from_aggregate <- sweep(
    as.matrix(aggregate_expression),
    MARGIN = 2L,
    STATS = domain_spot_counts,
    FUN = "/"
)
average_matrix <- as.matrix(average_expression)
maximum_absolute_difference <- max(
    abs(mean_from_aggregate - average_matrix),
    na.rm = TRUE
)
if (!is.finite(maximum_absolute_difference) ||
    maximum_absolute_difference > 1e-8) {
    stop(
        "AggregateExpression / n_spots does not match AverageExpression. ",
        "Maximum absolute difference: ", maximum_absolute_difference
    )
}

write.csv(
    data.frame(
        domain = domain_levels,
        n_spots = unname(domain_spot_counts)
    ),
    file.path(output_dir, "domain_spot_counts.csv"),
    row.names = FALSE,
    quote = FALSE
)
write_expression_matrix(
    aggregate_expression,
    file.path(
        output_dir,
        "RNA_AggregateExpression_sum_counts_by_domain.csv.gz"
    )
)
write_expression_matrix(
    average_expression,
    file.path(
        output_dir,
        "RNA_AverageExpression_mean_counts_by_domain.csv.gz"
    )
)
write.csv(
    data.frame(
        input_rds = basename(input_rds),
        assay = "RNA",
        layer = "counts",
        grouping_field = domain_column,
        retained_spots = ncol(domain_object),
        retained_genes = nrow(domain_object),
        maximum_absolute_difference = maximum_absolute_difference,
        stringsAsFactors = FALSE
    ),
    file.path(output_dir, "RNA_Aggregate_vs_Average_validation.csv"),
    row.names = FALSE,
    quote = TRUE
)

writeLines(
    capture.output(sessionInfo()),
    file.path(
        session_dir,
        "11_TO_GCN_domain_expression_export_sessionInfo.txt"
    )
)

message("Domain-level RNA expression export completed.")
message("Domains: ", paste(domain_levels, collapse = " -> "))
message("Retained spots: ", ncol(domain_object))
message("Genes: ", nrow(domain_object))
message("Maximum aggregate/average difference: ", maximum_absolute_difference)
message("Output directory: ", output_dir)
