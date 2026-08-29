# Plot QC and Harmony integration panels from the maize scRNA-seq reference
#
# Preferred lightweight input:
#   data/processed/sce_ref.rds
#
# Alternative input:
#   an in-memory Seurat object named sc_merged_filter_SCT2_inte, or
#   data/processed/sc_merged_filter_SCT2_inte.rds
#
# The script only reads the reference object. It does not modify or resave it.

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(patchwork)
})

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

first_existing <- function(candidates, available, description) {
  selected <- candidates[candidates %in% available]
  if (!length(selected)) {
    stop("Could not identify ", description, ". Available: ",
         paste(available, collapse = ", "))
  }
  selected[[1L]]
}

project_root <- find_project_root()
sce_rds <- file.path(project_root, "data", "processed", "sce_ref.rds")
seurat_rds <- file.path(
  project_root, "data", "processed", "sc_merged_filter_SCT2_inte.rds"
)
figure_dir <- file.path(
  project_root, "results", "figures", "10_scRNA_reference_integration"
)
table_dir <- file.path(
  project_root, "results", "tables", "10_scRNA_reference_integration"
)
session_dir <- file.path(project_root, "results", "sessionInfo")

dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(session_dir, recursive = TRUE, showWarnings = FALSE)

if (exists("sc_merged_filter_SCT2_inte", envir = .GlobalEnv, inherits = FALSE)) {
  reference <- get("sc_merged_filter_SCT2_inte", envir = .GlobalEnv)
  input_label <- "in-memory sc_merged_filter_SCT2_inte"
} else if (file.exists(sce_rds)) {
  reference <- readRDS(sce_rds)
  input_label <- "sce_ref.rds"
} else if (file.exists(seurat_rds)) {
  reference <- readRDS(seurat_rds)
  input_label <- "sc_merged_filter_SCT2_inte.rds"
} else {
  stop("No supported reference object was found.")
}

if (inherits(reference, "Seurat")) {
  metadata <- reference[[]]
  reduction_names <- SeuratObject::Reductions(reference)
  umap_name <- first_existing(
    c("umapharmony", "umap.harmony", "harmony.umap", "umap"),
    reduction_names,
    "Harmony UMAP reduction"
  )
  umap <- SeuratObject::Embeddings(reference, reduction = umap_name)
} else if (inherits(reference, "SingleCellExperiment")) {
  metadata <- as.data.frame(SummarizedExperiment::colData(reference))
  reduction_names <- SingleCellExperiment::reducedDimNames(reference)
  umap_name <- first_existing(
    c("UMAP.HARMONY", "umapharmony", "UMAP", "umap"),
    reduction_names,
    "Harmony UMAP reduction"
  )
  umap <- SingleCellExperiment::reducedDim(reference, umap_name)
} else {
  stop("Input must be a Seurat or SingleCellExperiment object.")
}

# Preserve the exact cell order when metadata and embeddings carry cell names.
if (nrow(metadata) != nrow(umap)) {
  stop(
    "Metadata and Harmony UMAP contain different numbers of cells: ",
    nrow(metadata), " versus ", nrow(umap), "."
  )
}
if (!is.null(rownames(metadata)) && !is.null(rownames(umap))) {
  if (!setequal(rownames(metadata), rownames(umap))) {
    stop("Metadata and Harmony UMAP do not contain the same cell identifiers.")
  }
  umap <- umap[rownames(metadata), , drop = FALSE]
}

metadata_columns <- colnames(metadata)
library_column <- first_existing(
  c("batch", "SRA_run", "library_id", "orig.ident", "sample_id"),
  metadata_columns,
  "library/batch column"
)
cluster_column <- first_existing(
  c("harmony_clusters", "seurat_clusters", "SCT_snn_res.2"),
  metadata_columns,
  "Harmony cluster column"
)

qc_features <- c("nFeature_RNA", "nCount_RNA")
mitochondrial_field <- intersect(
  c("percent.mt", "percent.mito"), metadata_columns
)[1L]
plastid_field <- intersect(
  c("percent.pl", "percent.pltd", "percent.chloroplast"), metadata_columns
)[1L]
qc_features <- unique(c(qc_features, mitochondrial_field, plastid_field))
qc_features <- qc_features[!is.na(qc_features) & qc_features %in% metadata_columns]
if (length(qc_features) < 2L) stop("Fewer than two QC fields were found.")

plot_data <- metadata
plot_data$plot_library <- factor(metadata[[library_column]])
plot_data$plot_cluster <- factor(metadata[[cluster_column]])
plot_data$UMAP_1 <- umap[, 1]
plot_data$UMAP_2 <- umap[, 2]

# Panel A: violin plots of genes, UMIs, mitochondrial reads, and plastid reads.
qc_long <- do.call(rbind, lapply(qc_features, function(feature) {
  data.frame(
    library = plot_data$plot_library,
    metric = feature,
    value = plot_data[[feature]],
    stringsAsFactors = FALSE
  )
}))
qc_long$metric <- factor(qc_long$metric, levels = qc_features)

panel_a <- ggplot(qc_long, aes(x = library, y = value)) +
  geom_violin(fill = "grey90", color = "grey25", linewidth = 0.3,
              scale = "width", trim = TRUE) +
  geom_jitter(width = 0.15, size = 0.10, alpha = 0.18, color = "black") +
  facet_wrap(~metric, nrow = 1, scales = "free_y") +
  labs(x = "Library", y = NULL) +
  theme_classic(base_size = 10) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold")
  )

ggsave(file.path(figure_dir, "scRNA_reference_QC_panel_A.png"), panel_a,
       width = max(8, 3.0 * length(qc_features)), height = 4.1, dpi = 300)
ggsave(file.path(figure_dir, "scRNA_reference_QC_panel_A.pdf"), panel_a,
       width = max(8, 3.0 * length(qc_features)), height = 4.1)

# Panel B: Harmony UMAP colored by library and by unsupervised cluster.
library_plot <- ggplot(
  plot_data, aes(UMAP_1, UMAP_2, color = plot_library)
) +
  geom_point(size = 0.25, alpha = 0.8) +
  labs(title = "Library", color = "Library",
       x = "Harmony UMAP-1", y = "Harmony UMAP-2") +
  guides(color = guide_legend(override.aes = list(size = 2, alpha = 1))) +
  theme_classic(base_size = 11)

cluster_centers <- plot_data |>
  dplyr::group_by(plot_cluster) |>
  dplyr::summarise(
    UMAP_1 = median(UMAP_1),
    UMAP_2 = median(UMAP_2),
    .groups = "drop"
  )

cluster_plot <- ggplot(
  plot_data, aes(UMAP_1, UMAP_2, color = plot_cluster)
) +
  geom_point(size = 0.25, alpha = 0.8) +
  geom_text(
    data = cluster_centers,
    aes(label = plot_cluster),
    color = "black", fontface = "bold", size = 3, show.legend = FALSE
  ) +
  labs(title = "Harmony clusters", color = "Cluster",
       x = "Harmony UMAP-1", y = "Harmony UMAP-2") +
  guides(color = guide_legend(override.aes = list(size = 2, alpha = 1), ncol = 2)) +
  theme_classic(base_size = 11)

panel_b <- library_plot | cluster_plot

ggsave(file.path(figure_dir, "scRNA_reference_Harmony_panel_B.png"), panel_b,
       width = 12, height = 5.7, dpi = 300)
ggsave(file.path(figure_dir, "scRNA_reference_Harmony_panel_B.pdf"), panel_b,
       width = 12, height = 5.7)

combined_figure <- (panel_a / panel_b) +
  plot_layout(heights = c(0.8, 1.2)) +
  plot_annotation(tag_levels = "A")

ggsave(file.path(figure_dir, "Figure_scRNA_reference_QC_and_Harmony.png"),
       combined_figure, width = 13, height = 10, dpi = 300)
ggsave(file.path(figure_dir, "Figure_scRNA_reference_QC_and_Harmony.pdf"),
       combined_figure, width = 13, height = 10)

qc_summary <- plot_data |>
  dplyr::group_by(plot_library) |>
  dplyr::summarise(
    n_cells = dplyr::n(),
    dplyr::across(
      dplyr::all_of(qc_features),
      list(
        median = ~median(.x, na.rm = TRUE),
        q05 = ~quantile(.x, 0.05, na.rm = TRUE),
        q95 = ~quantile(.x, 0.95, na.rm = TRUE)
      ),
      .names = "{.col}_{.fn}"
    ),
    .groups = "drop"
  )
write.csv(qc_summary,
          file.path(table_dir, "scRNA_reference_QC_summary_by_library.csv"),
          row.names = FALSE)

cluster_summary <- plot_data |>
  dplyr::count(plot_library, plot_cluster, name = "n_cells")
write.csv(cluster_summary,
          file.path(table_dir, "scRNA_reference_cells_by_library_and_cluster.csv"),
          row.names = FALSE)

plot_inputs <- data.frame(
  field = c("input", "library_column", "cluster_column", "UMAP_reduction", "QC_fields"),
  value = c(input_label, library_column, cluster_column, umap_name,
            paste(qc_features, collapse = ","))
)
write.csv(plot_inputs,
          file.path(table_dir, "scRNA_reference_plot_inputs.csv"),
          row.names = FALSE)

writeLines(capture.output(sessionInfo()),
           file.path(session_dir, "10_scRNA_reference_QC_plots_sessionInfo.txt"))

message("scRNA reference QC and Harmony figures completed from ", input_label, ".")
message("Panel A: ", file.path(figure_dir, "scRNA_reference_QC_panel_A.png"))
message("Panel B: ", file.path(figure_dir, "scRNA_reference_Harmony_panel_B.png"))
message("Combined: ", file.path(figure_dir, "Figure_scRNA_reference_QC_and_Harmony.png"))
