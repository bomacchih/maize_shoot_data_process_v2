# SCINA annotation and cell-type UMAP validation for the maize scRNA reference
#
# Marker input (CSV):
#   data/metadata/scRNA_reference/SCINA_marker_table.csv
# Required columns:
#   gene_id, cell_type
# Optional ranking columns (first available is used):
#   avg_log2FC, avg_logFC, marker_rank, rank
#
# Reference inputs, in order of preference:
#   1. in-memory Seurat object: sc_merged_filter_SCT2_inte
#   2. data/processed/sc_merged_filter_SCT2_inte.rds
#   3. data/processed/sce_ref.rds (plot-only when the full object is unavailable)
#
# A new RDS is written only when SCINA is run:
#   data/processed/sc_merged_filter_SCT2_inte_SCINA.rds
# The input object and its active identities are not modified.

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(patchwork)
})

set.seed(1)

find_project_root <- function(path = getwd()) {
  path <- normalizePath(path, winslash = "/", mustWork = TRUE)
  repeat {
    if (dir.exists(file.path(path, "data")) &&
        dir.exists(file.path(path, "scripts"))) return(path)
    parent <- dirname(path)
    if (identical(parent, path)) stop("Project root not found.")
    path <- parent
  }
}

first_existing <- function(candidates, available, description,
                           required = TRUE) {
  selected <- candidates[candidates %in% available]
  if (length(selected)) return(selected[[1L]])
  if (required) {
    stop("Could not identify ", description, ". Available: ",
         paste(available, collapse = ", "))
  }
  NA_character_
}

prepare_signatures <- function(marker_table, available_genes,
                               maximum_markers = 100L) {
  required_columns <- c("gene_id", "cell_type")
  missing_columns <- setdiff(required_columns, colnames(marker_table))
  if (length(missing_columns)) {
    stop("Marker table is missing: ", paste(missing_columns, collapse = ", "))
  }

  rank_column <- first_existing(
    c("avg_log2FC", "avg_logFC", "marker_rank", "rank"),
    colnames(marker_table), "marker-ranking column", required = FALSE
  )

  marker_table <- marker_table |>
    dplyr::mutate(
      gene_id = trimws(as.character(gene_id)),
      cell_type = trimws(as.character(cell_type)),
      input_order = seq_len(dplyr::n())
    ) |>
    dplyr::filter(!is.na(gene_id), nzchar(gene_id),
                  !is.na(cell_type), nzchar(cell_type),
                  gene_id %in% available_genes) |>
    dplyr::distinct(cell_type, gene_id, .keep_all = TRUE)

  if (!is.na(rank_column)) {
    decreasing_rank <- rank_column %in% c("avg_log2FC", "avg_logFC")
    rank_values <- suppressWarnings(as.numeric(marker_table[[rank_column]]))
    marker_table$rank_value <- if (decreasing_rank) -rank_values else rank_values
    marker_table <- marker_table |>
      dplyr::arrange(cell_type, is.na(rank_value), rank_value, input_order)
  } else {
    marker_table <- marker_table |>
      dplyr::arrange(cell_type, input_order)
  }

  # Remove markers assigned to more than one cell type.
  shared_markers <- marker_table |>
    dplyr::distinct(cell_type, gene_id) |>
    dplyr::count(gene_id, name = "n_cell_types") |>
    dplyr::filter(n_cell_types > 1L) |>
    dplyr::pull(gene_id)

  marker_table <- marker_table |>
    dplyr::filter(!gene_id %in% shared_markers) |>
    dplyr::group_by(cell_type) |>
    dplyr::slice_head(n = maximum_markers) |>
    dplyr::ungroup() |>
    dplyr::add_count(cell_type, name = "n_markers") |>
    dplyr::filter(n_markers >= 2L)

  signatures <- split(marker_table$gene_id, marker_table$cell_type)
  signatures <- lapply(signatures, unique)
  if (!length(signatures)) {
    stop("No cell type retained at least two non-overlapping marker genes.")
  }

  list(
    signatures = signatures,
    retained_table = marker_table,
    shared_markers = shared_markers,
    rank_column = rank_column
  )
}

standardize_scina_probabilities <- function(probabilities, cell_names) {
  probabilities <- as.matrix(probabilities)
  if (nrow(probabilities) == length(cell_names)) {
    rownames(probabilities) <- cell_names
    return(probabilities)
  }
  if (ncol(probabilities) == length(cell_names)) {
    probabilities <- t(probabilities)
    rownames(probabilities) <- cell_names
    return(probabilities)
  }
  stop("SCINA posterior matrix cannot be aligned to the input cells.")
}

project_root <- find_project_root()
marker_csv <- file.path(
  project_root, "data", "metadata", "scRNA_reference",
  "SCINA_marker_table.csv"
)
sce_rds <- file.path(project_root, "data", "processed", "sce_ref.rds")
seurat_rds <- file.path(
  project_root, "data", "processed", "sc_merged_filter_SCT2_inte.rds"
)
annotated_rds <- file.path(
  project_root, "data", "processed", "sc_merged_filter_SCT2_inte_SCINA.rds"
)
figure_dir <- file.path(
  project_root, "results", "figures", "10_scRNA_SCINA_annotation"
)
table_dir <- file.path(
  project_root, "results", "tables", "10_scRNA_SCINA_annotation"
)
session_dir <- file.path(project_root, "results", "sessionInfo")

dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(session_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(annotated_rds), recursive = TRUE, showWarnings = FALSE)

if (exists("sc_merged_filter_SCT2_inte", envir = .GlobalEnv,
           inherits = FALSE)) {
  reference <- get("sc_merged_filter_SCT2_inte", envir = .GlobalEnv)
  input_label <- "in-memory sc_merged_filter_SCT2_inte"
} else if (file.exists(seurat_rds)) {
  reference <- readRDS(seurat_rds)
  input_label <- "sc_merged_filter_SCT2_inte.rds"
} else if (file.exists(sce_rds)) {
  reference <- readRDS(sce_rds)
  input_label <- "sce_ref.rds (plot-only fallback)"
} else {
  stop("No supported reference object was found.")
}

is_seurat <- inherits(reference, "Seurat")
is_sce <- inherits(reference, "SingleCellExperiment")
if (!is_seurat && !is_sce) {
  stop("Reference must be a Seurat or SingleCellExperiment object.")
}

scina_was_run <- FALSE
annotation_column <- NA_character_

if (is_seurat && file.exists(marker_csv)) {
  if (!requireNamespace("SCINA", quietly = TRUE)) {
    stop("SCINA is required. Install SCINA v1.2.0 before running this step.")
  }
  if (utils::packageVersion("SCINA") != "1.2.0") {
    warning("Validated with SCINA 1.2.0; installed version is ",
            as.character(utils::packageVersion("SCINA")), ".")
  }
  if (!"SCT" %in% SeuratObject::Assays(reference)) {
    stop("The Seurat object does not contain an SCT assay.")
  }

  original_idents <- SeuratObject::Idents(reference)
  marker_table <- read.csv(marker_csv, check.names = FALSE,
                           stringsAsFactors = FALSE)
  available_genes <- rownames(reference[["SCT"]])
  prepared <- prepare_signatures(marker_table, available_genes, 100L)

  marker_genes <- unique(unlist(prepared$signatures, use.names = FALSE))
  expression_sparse <- SeuratObject::GetAssayData(
    reference, assay = "SCT", layer = "data"
  )[marker_genes, , drop = FALSE]
  expression_matrix <- as.matrix(expression_sparse)
  expression_matrix[!is.finite(expression_matrix)] <- 0

  set.seed(1)
  scina_result <- SCINA::SCINA(
    expression_matrix,
    prepared$signatures,
    max_iter = 100,
    convergence_n = 10,
    sensitivity_cutoff = 1,
    rm_overlap = TRUE,
    allow_unknown = TRUE
  )

  posterior <- standardize_scina_probabilities(
    scina_result$probabilities, colnames(expression_matrix)
  )
  posterior <- posterior[colnames(reference), , drop = FALSE]
  if (is.null(colnames(posterior))) {
    stop("SCINA posterior matrix does not contain cell-type column names.")
  }
  maximum_index <- max.col(posterior, ties.method = "first")
  maximum_posterior <- apply(posterior, 1, max, na.rm = TRUE)
  scina_labels <- colnames(posterior)[maximum_index]
  names(scina_labels) <- rownames(posterior)
  scina_labels[!is.finite(maximum_posterior) | maximum_posterior < 0.5] <- "Unknown"

  reference$celltype_scina <- factor(scina_labels)
  reference$celltype_scina_max_posterior <- maximum_posterior
  SeuratObject::Idents(reference) <- original_idents
  stopifnot(identical(SeuratObject::Idents(reference), original_idents))

  write.csv(
    prepared$retained_table,
    file.path(table_dir, "SCINA_marker_signatures_retained.csv"),
    row.names = FALSE
  )
  write.csv(
    data.frame(gene_id = prepared$shared_markers),
    file.path(table_dir, "SCINA_shared_markers_removed.csv"),
    row.names = FALSE
  )
  write.csv(
    data.frame(
      cell = colnames(reference),
      celltype_scina = scina_labels,
      maximum_posterior = maximum_posterior
    ),
    file.path(table_dir, "SCINA_cell_assignments.csv"),
    row.names = FALSE
  )

  saveRDS(reference, annotated_rds, compress = FALSE)
  annotation_column <- "celltype_scina"
  scina_was_run <- TRUE
} else {
  metadata_columns <- if (is_seurat) {
    colnames(reference[[]])
  } else {
    colnames(SummarizedExperiment::colData(reference))
  }
  annotation_column <- first_existing(
    c("celltype_scina", "celltype_scina_histo", "predicted_cell_type",
      "ident", "label"),
    metadata_columns,
    "existing cell-type annotation"
  )
  if (is_seurat && !file.exists(marker_csv)) {
    message("Marker CSV not found; plotting existing annotations from ",
            annotation_column, ".")
  }
}

# Extract annotation metadata and Harmony UMAP coordinates for Figure C.
if (is_seurat) {
  plot_metadata <- reference[[]]
  reductions <- SeuratObject::Reductions(reference)
  umap_name <- first_existing(
    c("umap_harmony", "umapharmony", "umap.harmony", "harmony.umap",
      "umap"),
    reductions, "Harmony UMAP reduction"
  )
  umap <- SeuratObject::Embeddings(reference, reduction = umap_name)
} else {
  plot_metadata <- as.data.frame(SummarizedExperiment::colData(reference))
  reductions <- SingleCellExperiment::reducedDimNames(reference)
  umap_name <- first_existing(
    c("UMAP_HARMONY", "UMAP.HARMONY", "umap_harmony", "umapharmony",
      "UMAP", "umap"),
    reductions, "Harmony UMAP reduction"
  )
  umap <- SingleCellExperiment::reducedDim(reference, umap_name)
}

plot_data <- data.frame(
  UMAP_1 = umap[, 1],
  UMAP_2 = umap[, 2],
  cell_type = as.character(plot_metadata[[annotation_column]]),
  stringsAsFactors = FALSE
)
plot_data$cell_type[is.na(plot_data$cell_type) |
                      !nzchar(plot_data$cell_type)] <- "Unknown"
cell_types <- sort(unique(plot_data$cell_type))

highlight_plots <- lapply(cell_types, function(current_type) {
  ggplot(plot_data, aes(UMAP_1, UMAP_2)) +
    geom_point(color = "grey82", size = 0.25, alpha = 0.7) +
    geom_point(
      data = plot_data[plot_data$cell_type == current_type, , drop = FALSE],
      color = "#E41A1C", size = 0.35, alpha = 0.9
    ) +
    labs(title = current_type, x = "Harmony UMAP-1", y = "Harmony UMAP-2") +
    coord_equal() +
    theme_classic(base_size = 9) +
    theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 9))
})

figure_c <- wrap_plots(highlight_plots, ncol = 4) +
  plot_annotation(
    title = "C   UMAP highlights of SCINA cell-type annotations"
  )

figure_height <- max(6, ceiling(length(cell_types) / 4) * 3.1)
ggsave(
  file.path(figure_dir, "scRNA_reference_SCINA_celltype_panel_C.png"),
  figure_c, width = 12, height = figure_height, dpi = 300
)
ggsave(
  file.path(figure_dir, "scRNA_reference_SCINA_celltype_panel_C.pdf"),
  figure_c, width = 12, height = figure_height
)

annotation_summary <- as.data.frame(sort(table(plot_data$cell_type),
                                         decreasing = TRUE))
colnames(annotation_summary) <- c("cell_type", "n_cells")
write.csv(
  annotation_summary,
  file.path(table_dir, "SCINA_celltype_summary.csv"),
  row.names = FALSE
)

run_record <- data.frame(
  field = c("input", "SCINA_was_run", "annotation_column",
            "UMAP_reduction", "seed", "posterior_cutoff"),
  value = c(input_label, scina_was_run, annotation_column,
            umap_name, 1, 0.5)
)
write.csv(run_record, file.path(table_dir, "SCINA_run_record.csv"),
          row.names = FALSE)

writeLines(
  capture.output(sessionInfo()),
  file.path(session_dir, "11_scRNA_SCINA_annotation_sessionInfo.txt")
)

message("SCINA annotation/validation workflow completed from ", input_label, ".")
message("Annotation field: ", annotation_column)
message("Figure C: ", file.path(
  figure_dir, "scRNA_reference_SCINA_celltype_panel_C.png"
))
if (scina_was_run) message("Annotated RDS: ", annotated_rds)
