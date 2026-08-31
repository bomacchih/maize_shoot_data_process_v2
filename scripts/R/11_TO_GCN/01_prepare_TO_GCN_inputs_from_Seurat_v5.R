#!/usr/bin/env Rscript

# Prepare single-time-series TO-GCN inputs from the embryonic-leaf Seurat v5
# object. Raw RNA UMIs are summed within each structural domain and divided by
# the number of retained spots, yielding mean UMI per spot.
#
# Run from the repository root:
#   source("scripts/R/11_TO_GCN/01_prepare_TO_GCN_inputs_from_Seurat_v5.R")
#
# Main input:
#   data/processed/XGE202122_S5_subset_embleaf_harmony_join.rds
#
# Required TF annotation (use either format):
#   data/reference/TO_GCN/maize_TF_annotation.csv
#     Required columns: gene_id, gene_name, tf_family
#   data/reference/TO_GCN/gene_expressions_v3.xlsx
#     Recovered workbook; sheet TF_genes supplies the same annotation and can
#     also be used to validate the archived L1-L13 assignments.
#
# Main outputs:
#   results/tables/11_TO_GCN/input/domain_spot_counts.csv
#   results/tables/11_TO_GCN/input/all_genes_mean_UMI_per_spot.csv.gz
#   results/tables/11_TO_GCN/input/TF_genes_mean_UMI_per_spot.csv
#   results/tables/11_TO_GCN/input/TF_genes_gene_wise_zscore.csv
#   results/tables/11_TO_GCN/input/All_gene_expression.tsv
#   results/tables/11_TO_GCN/input/TF_expression.tsv
#   results/tables/11_TO_GCN/input/seeds.txt
#   results/sessionInfo/11_TO_GCN_prepare_sessionInfo.txt

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

write_matrix_csv <- function(matrix_to_write, path, compress = FALSE) {
    output <- data.frame(
        gene_id = rownames(matrix_to_write),
        as.data.frame(matrix_to_write, check.names = FALSE),
        check.names = FALSE
    )
    connection <- if (compress) gzfile(path, open = "wt") else path
    if (compress) on.exit(close(connection), add = TRUE)
    write.csv(output, connection, row.names = FALSE, quote = FALSE)
}

write_togcn_tsv <- function(matrix_to_write, path) {
    output <- data.frame(
        gene_id = rownames(matrix_to_write),
        as.data.frame(matrix_to_write, check.names = FALSE),
        check.names = FALSE
    )
    # The original C++ programs parse whitespace-delimited records with no
    # header. Keep full numeric precision and prohibit missing values.
    if (anyNA(output)) stop("TO-GCN input contains missing values: ", path)
    write.table(
        output,
        file = path,
        sep = "\t",
        row.names = FALSE,
        col.names = FALSE,
        quote = FALSE,
        na = ""
    )
}

read_tf_annotation <- function(csv_path, workbook_path) {
    recovered_levels <- NULL

    if (file.exists(csv_path)) {
        annotation <- read.csv(
            csv_path,
            stringsAsFactors = FALSE,
            check.names = FALSE
        )
    } else if (file.exists(workbook_path)) {
        if (!requireNamespace("readxl", quietly = TRUE)) {
            stop(
                "The recovered TF workbook was found, but package 'readxl' is absent. ",
                "Install readxl or provide maize_TF_annotation.csv."
            )
        }
        annotation <- as.data.frame(
            readxl::read_excel(workbook_path, sheet = "TF_genes"),
            stringsAsFactors = FALSE,
            check.names = FALSE
        )
    } else {
        stop(
            "No TF annotation was found. Provide either:\n  ", csv_path,
            "\nor:\n  ", workbook_path
        )
    }

    gene_column <- detect_column(
        names(annotation),
        c("gene_id", "gene", "tf_gene_id")
    )
    name_column <- detect_column(
        names(annotation),
        c("gene_name", "genename", "symbol")
    )
    family_column <- detect_column(
        names(annotation),
        c("tf_family", "gene_family", "family")
    )
    level_column <- detect_column(
        names(annotation),
        c("archived_level", "level_in_gcn", "cluster_level", "level", "cluster")
    )

    if (is.na(gene_column)) {
        stop("The TF annotation does not contain a recognizable gene-ID column.")
    }

    output <- data.frame(
        gene_id = trimws(as.character(annotation[[gene_column]])),
        gene_name = if (is.na(name_column)) NA_character_ else
            trimws(as.character(annotation[[name_column]])),
        tf_family = if (is.na(family_column)) NA_character_ else
            trimws(as.character(annotation[[family_column]])),
        stringsAsFactors = FALSE
    )
    output <- output[!is.na(output$gene_id) & nzchar(output$gene_id), , drop = FALSE]
    output <- output[!duplicated(output$gene_id), , drop = FALSE]

    if (!is.na(level_column)) {
        recovered_levels <- data.frame(
            gene_id = trimws(as.character(annotation[[gene_column]])),
            recovered_level = suppressWarnings(as.integer(annotation[[level_column]])),
            stringsAsFactors = FALSE
        )
        recovered_levels <- recovered_levels[
            !is.na(recovered_levels$gene_id) &
                nzchar(recovered_levels$gene_id) &
                !is.na(recovered_levels$recovered_level),
            ,
            drop = FALSE
        ]
        recovered_levels <- recovered_levels[!duplicated(recovered_levels$gene_id), ]
    }

    list(annotation = output, recovered_levels = recovered_levels)
}

get_raw_rna_counts <- function(object) {
    if (!"RNA" %in% Assays(object)) {
        stop("The input object does not contain an RNA assay.")
    }

    layer_names <- Layers(object[["RNA"]])
    count_layers <- grep("^counts($|\\.)", layer_names, value = TRUE)
    if (length(count_layers) == 0L) {
        stop("No raw RNA count layer was found. Available layers: ",
             paste(layer_names, collapse = ", "))
    }

    matrices <- lapply(count_layers, function(layer_name) {
        x <- LayerData(object, assay = "RNA", layer = layer_name)
        if (!inherits(x, "sparseMatrix")) x <- as(x, "dgCMatrix")
        x
    })

    if (length(matrices) == 1L) return(matrices[[1L]])

    common_features <- Reduce(intersect, lapply(matrices, rownames))
    if (length(common_features) == 0L) {
        stop("The RNA count layers have no common features.")
    }
    matrices <- lapply(
        matrices,
        function(x) x[common_features, , drop = FALSE]
    )
    combined <- do.call(cbind, matrices)
    if (anyDuplicated(colnames(combined))) {
        stop("Duplicated spot names were found while combining RNA count layers.")
    }
    combined
}

project_root <- find_project_root()

input_rds <- file.path(
    project_root,
    "data", "processed", "XGE202122_S5_subset_embleaf_harmony_join.rds"
)
reference_dir <- file.path(project_root, "data", "reference", "TO_GCN")
tf_annotation_csv <- file.path(reference_dir, "maize_TF_annotation.csv")
recovered_workbook <- file.path(reference_dir, "gene_expressions_v3.xlsx")
output_dir <- file.path(project_root, "results", "tables", "11_TO_GCN", "input")
session_dir <- file.path(project_root, "results", "sessionInfo")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(session_dir, recursive = TRUE, showWarnings = FALSE)

domain_levels <- c("SAM", "P1_P2", "P3", "P4", "P5")
domain_column_preference <- c("domains", "domain", "structural_domain")

# Recovered candidate from the earlier TO-GCN workbook. Change this vector if
# the original seed record or a new biological hypothesis specifies otherwise.
seed_gene_ids <- c("Zm00001eb222520") # mybr90
positive_pcc_cutoff <- 0.95

if (!file.exists(input_rds)) stop("Input RDS not found: ", input_rds)
embryonic_leaf <- readRDS(input_rds)
raw_counts <- get_raw_rna_counts(embryonic_leaf)

metadata <- embryonic_leaf[[]]
domain_column <- detect_column(names(metadata), domain_column_preference)
if (is.na(domain_column)) {
    stop(
        "Could not identify structural-domain metadata. Available columns: ",
        paste(names(metadata), collapse = ", ")
    )
}

missing_metadata_spots <- setdiff(colnames(raw_counts), rownames(metadata))
if (length(missing_metadata_spots) > 0L) {
    stop(length(missing_metadata_spots), " raw-count spots are missing from metadata.")
}
metadata <- metadata[colnames(raw_counts), , drop = FALSE]
domain_values <- trimws(as.character(metadata[[domain_column]]))
keep_spots <- !is.na(domain_values) & domain_values %in% domain_levels

raw_counts <- raw_counts[, keep_spots, drop = FALSE]
domain_values <- factor(domain_values[keep_spots], levels = domain_levels)
if (ncol(raw_counts) == 0L) stop("No spots remained in SAM through P5.")

domain_spot_counts <- as.integer(table(domain_values)[domain_levels])
if (any(domain_spot_counts == 0L)) {
    stop(
        "At least one required domain has no spots: ",
        paste(domain_levels[domain_spot_counts == 0L], collapse = ", ")
    )
}

membership <- sparseMatrix(
    i = seq_along(domain_values),
    j = as.integer(domain_values),
    x = 1,
    dims = c(length(domain_values), length(domain_levels)),
    dimnames = list(colnames(raw_counts), domain_levels)
)

domain_umi_sums <- raw_counts %*% membership
domain_mean_umi <- sweep(
    as.matrix(domain_umi_sums),
    MARGIN = 2L,
    STATS = domain_spot_counts,
    FUN = "/"
)
colnames(domain_mean_umi) <- domain_levels

total_umi <- Matrix::rowSums(raw_counts)
expressed <- is.finite(total_umi) & total_umi > 0
domain_mean_umi <- domain_mean_umi[expressed, , drop = FALSE]

# Constant profiles have undefined Pearson correlations and cannot contribute
# edges to TO-GCN. Retain them in the background table but exclude them from
# the C++ network input.
profile_sd <- apply(domain_mean_umi, 1L, sd)
network_eligible <- is.finite(profile_sd) & profile_sd > 0
network_all_genes <- domain_mean_umi[network_eligible, , drop = FALSE]

tf_info <- read_tf_annotation(tf_annotation_csv, recovered_workbook)
tf_annotation <- tf_info$annotation
tf_gene_ids <- intersect(tf_annotation$gene_id, rownames(network_all_genes))
if (length(tf_gene_ids) < 2L) {
    stop("Fewer than two annotated TF genes overlap the expressed-gene matrix.")
}

tf_annotation <- tf_annotation[match(tf_gene_ids, tf_annotation$gene_id), , drop = FALSE]
tf_mean_umi <- network_all_genes[tf_gene_ids, , drop = FALSE]

tf_zscore <- t(scale(t(tf_mean_umi), center = TRUE, scale = TRUE))
tf_zscore <- tf_zscore[apply(tf_zscore, 1L, function(x) all(is.finite(x))), , drop = FALSE]

missing_seeds <- setdiff(seed_gene_ids, rownames(tf_mean_umi))
if (length(missing_seeds) > 0L) {
    stop(
        "Seed gene(s) are absent from the TF input: ",
        paste(missing_seeds, collapse = ", ")
    )
}

write.csv(
    data.frame(domain = domain_levels, n_spots = domain_spot_counts),
    file.path(output_dir, "domain_spot_counts.csv"),
    row.names = FALSE,
    quote = FALSE
)
write_matrix_csv(
    domain_mean_umi,
    file.path(output_dir, "all_genes_mean_UMI_per_spot.csv.gz"),
    compress = TRUE
)
write_matrix_csv(
    tf_mean_umi,
    file.path(output_dir, "TF_genes_mean_UMI_per_spot.csv")
)
write_matrix_csv(
    tf_zscore,
    file.path(output_dir, "TF_genes_gene_wise_zscore.csv")
)
write.csv(
    tf_annotation,
    file.path(output_dir, "TF_annotation_used.csv"),
    row.names = FALSE,
    quote = FALSE
)
write_togcn_tsv(
    network_all_genes,
    file.path(output_dir, "All_gene_expression.tsv")
)
write_togcn_tsv(
    tf_mean_umi,
    file.path(output_dir, "TF_expression.tsv")
)
writeLines(seed_gene_ids, file.path(output_dir, "seeds.txt"))

write.csv(
    data.frame(
        parameter = c(
            "input_rds", "domain_metadata_column", "ordered_domains",
            "normalization", "number_of_all_expressed_genes",
            "number_of_network_eligible_genes", "number_of_TFs",
            "seed_gene_ids", "positive_PCC_cutoff"
        ),
        value = c(
            basename(input_rds), domain_column, paste(domain_levels, collapse = ";"),
            "sum raw RNA UMIs per domain / number of spots in domain",
            nrow(domain_mean_umi), nrow(network_all_genes), nrow(tf_mean_umi),
            paste(seed_gene_ids, collapse = ";"), positive_pcc_cutoff
        ),
        stringsAsFactors = FALSE
    ),
    file.path(output_dir, "TO_GCN_input_parameters.csv"),
    row.names = FALSE,
    quote = TRUE
)

if (!is.null(tf_info$recovered_levels)) {
    recovered_levels <- tf_info$recovered_levels[
        tf_info$recovered_levels$gene_id %in% tf_gene_ids,
        ,
        drop = FALSE
    ]
    write.csv(
        recovered_levels,
        file.path(output_dir, "recovered_reference_TF_levels.csv"),
        row.names = FALSE,
        quote = FALSE
    )
}

writeLines(
    capture.output(sessionInfo()),
    file.path(session_dir, "11_TO_GCN_prepare_sessionInfo.txt")
)

message("TO-GCN input preparation completed.")
message("Domains: ", paste(domain_levels, collapse = " -> "))
message("Expressed genes: ", nrow(domain_mean_umi))
message("Network-eligible genes: ", nrow(network_all_genes))
message("Annotated TFs: ", nrow(tf_mean_umi))
message("Output directory: ", output_dir)
