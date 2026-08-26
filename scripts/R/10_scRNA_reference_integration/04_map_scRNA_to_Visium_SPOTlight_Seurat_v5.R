# Map annotated scRNA-seq cell types onto Visium spots with Seurat and SPOTlight
#
# Main inputs:
#   data/processed/sc_merged_filter_SCT2_inte_SCINA.rds
#   data/processed/maize_shoot_14samples_SCT_harmony_seurat_v5.rds
#
# Recommended for large objects:
#   Load both objects in RStudio as `sc_reference` and `visium_query` before
#   sourcing this script. In-memory objects take precedence over RDS files.
#
# Output:
#   data/processed/maize_shoot_14samples_celltype_mapped_SPOTlight_seurat_v5.rds

suppressPackageStartupMessages({
  library(Seurat)
  library(SingleCellExperiment)
  library(SpatialExperiment)
  library(scuttle)
  library(scran)
  library(SPOTlight)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
})

set.seed(1)

# -----------------------------
# User-adjustable settings
# -----------------------------

run_hard_label_transfer <- TRUE
run_spotlight_deconvolution <- TRUE
sections_to_run <- NULL       # NULL runs every section; use c("VR03_S2") to test
reference_cells_per_type <- 100L
markers_per_type <- 100L
n_hvg <- 3000L
n_pcs <- 30L
spotlight_min_prop <- 0.01
high_purity_cutoff <- 0.60
vascular_display_cutoff <- 0.10
sam_display_cutoff <- 0.05

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

safe_spotlight_name <- function(cell_type) {
  paste0("SPOT_", make.names(cell_type))
}

standardize_cell_type <- function(x) {
  gsub("^SPOT_", "", as.character(x))
}

get_all_tissue_coordinates <- function(object) {
  image_names <- Images(object)
  if (!length(image_names)) stop("The Visium object contains no spatial images.")

  coordinate_list <- lapply(image_names, function(image_name) {
    coordinates <- GetTissueCoordinates(object, image = image_name)
    if (!nrow(coordinates)) return(NULL)
    coordinates$spatial_image <- image_name
    coordinates$spot_barcode <- rownames(coordinates)
    coordinates
  })
  coordinates <- bind_rows(coordinate_list)
  if (!nrow(coordinates)) stop("No tissue coordinates were recovered.")
  coordinates <- coordinates[!duplicated(coordinates$spot_barcode), , drop = FALSE]
  rownames(coordinates) <- coordinates$spot_barcode
  coordinates
}

select_xy_columns <- function(coordinates) {
  available <- colnames(coordinates)
  candidate_pairs <- list(
    c("imagecol", "imagerow"),
    c("x", "y"),
    c("col", "row"),
    c("pxl_col_in_fullres", "pxl_row_in_fullres")
  )
  for (pair in candidate_pairs) {
    if (all(pair %in% available)) return(pair)
  }
  numeric_columns <- available[vapply(coordinates, is.numeric, logical(1))]
  if (length(numeric_columns) < 2L) {
    stop("Could not identify two numeric spatial-coordinate columns.")
  }
  numeric_columns[1:2]
}

remove_organelle_and_ribosomal_genes <- function(
    genes, mitochondrial_genes = character(), plastid_genes = character()) {
  ribosomal_pattern <- "^(RPL|RPS|Rp[LlSs]|rpl|rps)"
  mitochondrial_pattern <- "^(MT-|mt-|mitochond)"
  plastid_pattern <- "^(PLTD-|pt-|chloroplast)"
  genes[
    !grepl(ribosomal_pattern, genes) &
      !grepl(mitochondrial_pattern, genes, ignore.case = TRUE) &
      !grepl(plastid_pattern, genes, ignore.case = TRUE) &
      !genes %in% mitochondrial_genes &
      !genes %in% plastid_genes
  ]
}

read_optional_gene_list <- function(path) {
  if (!file.exists(path)) return(character())
  unique(trimws(readLines(path, warn = FALSE))) |>
    (\(x) x[nzchar(x)])()
}

build_marker_table <- function(sce, genes_to_test, maximum_per_type = 100L) {
  marker_statistics <- scran::scoreMarkers(sce, subset.row = genes_to_test)
  marker_tables <- lapply(names(marker_statistics), function(cell_type) {
    marker_table <- as.data.frame(marker_statistics[[cell_type]])
    weight_column <- first_existing(
      c("mean.AUC", "median.AUC", "mean.logFC", "summary.logFC", "cohen"),
      colnames(marker_table), "marker effect-size/AUC column"
    )
    marker_table$gene <- rownames(marker_table)
    marker_table$cluster <- cell_type
    marker_table$weight <- suppressWarnings(as.numeric(marker_table[[weight_column]]))
    marker_table |>
      filter(is.finite(weight), gene %in% genes_to_test) |>
      arrange(desc(weight)) |>
      slice_head(n = maximum_per_type) |>
      select(gene, cluster, weight)
  })
  marker_table <- bind_rows(marker_tables)
  if (!nrow(marker_table)) stop("No SPOTlight marker genes were retained.")
  marker_table
}

downsample_reference <- function(sce, groups, maximum_per_type = 100L) {
  set.seed(1)
  indices <- split(seq_len(ncol(sce)), groups)
  retained <- unlist(lapply(indices, function(index) {
    sample(index, min(length(index), maximum_per_type))
  }), use.names = FALSE)
  sce[, retained, drop = FALSE]
}

normalize_rows <- function(matrix_object) {
  matrix_object <- as.matrix(matrix_object)
  matrix_object[!is.finite(matrix_object) | matrix_object < 0] <- 0
  row_totals <- rowSums(matrix_object)
  valid <- row_totals > 0
  matrix_object[valid, ] <- matrix_object[valid, , drop = FALSE] /
    row_totals[valid]
  matrix_object
}

project_root <- find_project_root()
reference_rds <- file.path(
  project_root, "data", "processed",
  "sc_merged_filter_SCT2_inte_SCINA.rds"
)
visium_rds <- file.path(
  project_root, "data", "processed",
  "maize_shoot_14samples_SCT_harmony_seurat_v5.rds"
)
output_rds <- file.path(
  project_root, "data", "processed",
  "maize_shoot_14samples_celltype_mapped_SPOTlight_seurat_v5.rds"
)
optional_module_marker_csv <- file.path(
  project_root, "data", "metadata", "scRNA_reference",
  "independent_celltype_module_markers.csv"
)
mitochondrial_gene_file <- file.path(
  project_root, "data", "metadata", "maize_mitochondrial_genes.txt"
)
plastid_gene_file <- file.path(
  project_root, "data", "metadata", "maize_plastid_genes.txt"
)
figure_dir <- file.path(
  project_root, "results", "figures", "12_scRNA_Visium_mapping"
)
table_dir <- file.path(
  project_root, "results", "tables", "12_scRNA_Visium_mapping"
)
session_dir <- file.path(project_root, "results", "sessionInfo")
section_pie_dir <- file.path(figure_dir, "section_scatterpies")

for (directory in c(dirname(output_rds), figure_dir, table_dir,
                    session_dir, section_pie_dir)) {
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
}

if (exists("sc_reference", envir = .GlobalEnv, inherits = FALSE)) {
  sc_reference <- get("sc_reference", envir = .GlobalEnv)
} else if (file.exists(reference_rds)) {
  sc_reference <- readRDS(reference_rds)
} else if (exists("sc_merged_filter_SCT2_inte", envir = .GlobalEnv,
                  inherits = FALSE)) {
  sc_reference <- get("sc_merged_filter_SCT2_inte", envir = .GlobalEnv)
} else {
  stop("Load an annotated scRNA reference as `sc_reference`, or create: ",
       reference_rds)
}

if (exists("visium_query", envir = .GlobalEnv, inherits = FALSE)) {
  visium_query <- get("visium_query", envir = .GlobalEnv)
} else if (file.exists(visium_rds)) {
  visium_query <- readRDS(visium_rds)
} else {
  stop("Load the combined Visium object as `visium_query`, or create: ",
       visium_rds)
}

stopifnot(inherits(sc_reference, "Seurat"), inherits(visium_query, "Seurat"))
original_visium_idents <- Idents(visium_query)
reference_metadata <- sc_reference[[]]
visium_metadata <- visium_query[[]]

celltype_column <- first_existing(
  c("celltype_scina", "celltype_scina_histo", "predicted_cell_type",
    "ident", "label"),
  colnames(reference_metadata), "scRNA cell-type annotation column"
)
reference_labels <- as.character(reference_metadata[[celltype_column]])
valid_reference_cells <- !is.na(reference_labels) & nzchar(reference_labels) &
  reference_labels != "Unknown"
sc_reference <- subset(sc_reference, cells = colnames(sc_reference)[valid_reference_cells])
sc_reference$celltype_mapping_reference <- factor(
  reference_labels[valid_reference_cells]
)

reference_assay <- first_existing(
  c("SCT", "RNA"), Assays(sc_reference), "reference expression assay"
)
query_assay <- first_existing(
  c("SCT", "RNA", "Spatial"), Assays(visium_query),
  "Visium expression assay"
)

# -----------------------------
# Hard label transfer with Seurat anchors and MapQuery
# -----------------------------

if (run_hard_label_transfer) {
  if (!"SCT" %in% Assays(sc_reference) || !"SCT" %in% Assays(visium_query)) {
    stop("SCT assays are required for the configured SCT anchor transfer.")
  }
  reference_reduction <- first_existing(
    c("integrated.harmony", "harmony", "pca"),
    Reductions(sc_reference), "reference PCA/Harmony reduction"
  )
  available_dimensions <- min(
    n_pcs,
    ncol(Embeddings(sc_reference, reduction = reference_reduction))
  )
  transfer_dimensions <- seq_len(available_dimensions)

  # MapQuery requires a UMAP model. Create a dedicated model if necessary.
  mapping_umap <- "reference.mapping.umap"
  sc_reference <- RunUMAP(
    sc_reference,
    reduction = reference_reduction,
    dims = transfer_dimensions,
    return.model = TRUE,
    reduction.name = mapping_umap,
    reduction.key = "refMAPUMAP_",
    verbose = FALSE
  )

  transfer_anchors <- FindTransferAnchors(
    reference = sc_reference,
    query = visium_query,
    normalization.method = "SCT",
    reference.assay = "SCT",
    query.assay = "SCT",
    reference.reduction = reference_reduction,
    dims = transfer_dimensions
  )

  visium_query <- MapQuery(
    anchorset = transfer_anchors,
    query = visium_query,
    reference = sc_reference,
    refdata = list(celltype = "celltype_mapping_reference"),
    reference.reduction = reference_reduction,
    reference.dims = transfer_dimensions,
    reduction.model = mapping_umap
  )

  # MapQuery creates `predicted.celltype` and `predicted.celltype.score`.
  if (!"predicted.celltype" %in% colnames(visium_query[[]])) {
    stop("MapQuery did not create predicted.celltype metadata.")
  }
}

# -----------------------------
# SPOTlight markers, HVGs, and reference downsampling
# -----------------------------

if (run_spotlight_deconvolution) {
  reference_counts <- GetAssayData(
    sc_reference, assay = reference_assay, layer = "counts"
  )
  sce_reference <- SingleCellExperiment(
    assays = list(counts = reference_counts),
    colData = S4Vectors::DataFrame(
      cell_type = sc_reference$celltype_mapping_reference
    )
  )
  sce_reference <- scuttle::logNormCounts(sce_reference)
  SingleCellExperiment::colLabels(sce_reference) <-
    factor(SummarizedExperiment::colData(sce_reference)$cell_type)

  mitochondrial_genes <- read_optional_gene_list(mitochondrial_gene_file)
  plastid_genes <- read_optional_gene_list(plastid_gene_file)
  eligible_genes <- remove_organelle_and_ribosomal_genes(
    rownames(sce_reference), mitochondrial_genes, plastid_genes
  )

  variance_model <- scran::modelGeneVar(
    sce_reference, subset.row = eligible_genes
  )
  hvg <- scran::getTopHVGs(variance_model, n = min(n_hvg, nrow(variance_model)))
  marker_table <- build_marker_table(
    sce_reference, eligible_genes, markers_per_type
  )
  write.csv(marker_table,
            file.path(table_dir, "SPOTlight_scRNA_marker_table.csv"),
            row.names = FALSE)
  write.csv(data.frame(gene = hvg),
            file.path(table_dir, "SPOTlight_scRNA_HVGs.csv"),
            row.names = FALSE)

  sce_reference_small <- downsample_reference(
    sce_reference,
    groups = SummarizedExperiment::colData(sce_reference)$cell_type,
    maximum_per_type = reference_cells_per_type
  )
  reference_cell_types <- levels(droplevels(
    factor(SummarizedExperiment::colData(sce_reference_small)$cell_type)
  ))

  spatial_assay <- first_existing(
    c("RNA", "Spatial"), Assays(visium_query), "raw Visium count assay"
  )
  spatial_counts <- GetAssayData(
    visium_query, assay = spatial_assay, layer = "counts"
  )
  all_coordinates <- get_all_tissue_coordinates(visium_query)
  xy_columns <- select_xy_columns(all_coordinates)

  section_column <- first_existing(
    c("section_id", "domain_section", "section"),
    colnames(visium_query[[]]), "section metadata column"
  )
  section_labels <- as.character(visium_query[[section_column, drop = TRUE]])
  names(section_labels) <- colnames(visium_query)
  available_sections <- sort(unique(section_labels[!is.na(section_labels) &
                                                      nzchar(section_labels)]))
  target_sections <- if (is.null(sections_to_run)) {
    available_sections
  } else {
    intersect(sections_to_run, available_sections)
  }
  if (!length(target_sections)) stop("No requested Visium sections were found.")

  all_proportions <- matrix(
    NA_real_, nrow = ncol(visium_query), ncol = length(reference_cell_types),
    dimnames = list(colnames(visium_query), reference_cell_types)
  )
  section_diagnostics <- vector("list", length(target_sections))
  names(section_diagnostics) <- target_sections

  for (section_id in target_sections) {
    message("Running SPOTlight for section: ", section_id)
    section_spots <- names(section_labels)[section_labels == section_id]
    section_spots <- intersect(
      section_spots,
      intersect(colnames(spatial_counts), rownames(all_coordinates))
    )
    if (length(section_spots) < 10L) {
      warning("Skipping ", section_id, ": fewer than 10 aligned tissue spots.")
      next
    }

    common_genes <- intersect(
      rownames(sce_reference_small), rownames(spatial_counts)
    )
    section_markers <- marker_table |>
      filter(gene %in% common_genes)
    section_hvg <- intersect(hvg, common_genes)
    if (!nrow(section_markers) || length(section_hvg) < 50L) {
      warning("Skipping ", section_id, ": insufficient shared markers/HVGs.")
      next
    }

    coordinates <- all_coordinates[section_spots, xy_columns, drop = FALSE]
    spatial_experiment <- SpatialExperiment(
      assays = list(
        counts = spatial_counts[common_genes, section_spots, drop = FALSE]
      ),
      colData = S4Vectors::DataFrame(
        visium_query[[]][section_spots, , drop = FALSE]
      ),
      spatialCoords = as.matrix(coordinates)
    )
    section_reference <- sce_reference_small[common_genes, , drop = FALSE]

    set.seed(1)
    spotlight_result <- SPOTlight::SPOTlight(
      x = section_reference,
      y = spatial_experiment,
      groups = as.character(
        SummarizedExperiment::colData(section_reference)$cell_type
      ),
      mgs = section_markers,
      hvg = section_hvg,
      weight_id = "weight",
      group_id = "cluster",
      gene_id = "gene",
      min_prop = spotlight_min_prop,
      slot_sc = "counts",
      slot_sp = "counts"
    )

    raw_proportions <- spotlight_result$mat
    cell_type_columns <- intersect(reference_cell_types,
                                   colnames(raw_proportions))
    proportions <- normalize_rows(
      raw_proportions[, cell_type_columns, drop = FALSE]
    )
    all_proportions[rownames(proportions), cell_type_columns] <- proportions

    section_diagnostics[[section_id]] <- data.frame(
      section_id = section_id,
      n_spots = nrow(proportions),
      n_shared_genes = length(common_genes),
      n_HVGs = length(section_hvg),
      n_markers = nrow(section_markers)
    )

    scatterpie_plot <- SPOTlight::plotSpatialScatterpie(
      # SPOTlight 1.16.0 requires a coordinate matrix for image-free plots.
      x = as.matrix(coordinates),
      y = proportions,
      cell_types = colnames(proportions),
      img = FALSE,
      scatterpie_alpha = 0.85,
      pie_scale = 0.45
    ) + ggtitle(paste("SPOTlight cell-type proportions:", section_id))
    ggsave(
      file.path(section_pie_dir, paste0(make.names(section_id), ".png")),
      scatterpie_plot, width = 8, height = 7, dpi = 300
    )
  }

  completed_spots <- rowSums(all_proportions, na.rm = TRUE) > 0
  proportion_matrix <- all_proportions[completed_spots, , drop = FALSE]
  if (!nrow(proportion_matrix)) {
    stop("SPOTlight did not return a positive proportion vector for any spot.")
  }
  proportion_matrix <- normalize_rows(proportion_matrix)
  all_proportions[completed_spots, ] <- proportion_matrix
  proportion_metadata <- all_proportions
  colnames(proportion_metadata) <- safe_spotlight_name(
    colnames(proportion_metadata)
  )
  visium_query <- AddMetaData(visium_query, proportion_metadata)

  maximum_index <- max.col(proportion_matrix, ties.method = "first")
  top_type <- colnames(proportion_matrix)[maximum_index]
  top_proportion <- apply(proportion_matrix, 1, max, na.rm = TRUE)
  entropy <- -rowSums(proportion_matrix * log(proportion_matrix + 1e-12))
  normalized_entropy <- if (ncol(proportion_matrix) > 1L) {
    entropy / log(ncol(proportion_matrix))
  } else {
    rep(0, length(entropy))
  }

  summary_metadata <- data.frame(
    SPOT_top_type = top_type,
    SPOT_top_prop = top_proportion,
    SPOT_high_purity = top_proportion >= high_purity_cutoff,
    SPOT_entropy = entropy,
    SPOT_entropy_normalized = normalized_entropy,
    row.names = rownames(proportion_matrix)
  )
  visium_query <- AddMetaData(visium_query, summary_metadata)

  write.csv(
    cbind(spot = rownames(proportion_matrix), as.data.frame(proportion_matrix)),
    file.path(table_dir, "SPOTlight_spot_celltype_proportions.csv"),
    row.names = FALSE
  )
  write.csv(
    bind_rows(section_diagnostics),
    file.path(table_dir, "SPOTlight_section_diagnostics.csv"),
    row.names = FALSE
  )
}

# -----------------------------
# Reliability diagnostics
# -----------------------------

metadata_final <- visium_query[[]]
if (all(c("predicted.celltype", "SPOT_top_type") %in%
        colnames(metadata_final))) {
  evaluable <- !is.na(metadata_final$predicted.celltype) &
    !is.na(metadata_final$SPOT_top_type)
  agreement_table <- data.frame(
    spot = rownames(metadata_final)[evaluable],
    transferred_label = as.character(metadata_final$predicted.celltype[evaluable]),
    SPOTlight_argmax = as.character(metadata_final$SPOT_top_type[evaluable]),
    agreement = as.character(metadata_final$predicted.celltype[evaluable]) ==
      as.character(metadata_final$SPOT_top_type[evaluable])
  )
  write.csv(agreement_table,
            file.path(table_dir, "SPOTlight_vs_Seurat_label_agreement.csv"),
            row.names = FALSE)
  write.csv(
    data.frame(
      n_evaluable_spots = nrow(agreement_table),
      agreement_rate = mean(agreement_table$agreement)
    ),
    file.path(table_dir, "SPOTlight_vs_Seurat_agreement_summary.csv"),
    row.names = FALSE
  )
}

# Optional independent module-score validation.
if (file.exists(optional_module_marker_csv)) {
  module_markers <- read.csv(optional_module_marker_csv,
                             stringsAsFactors = FALSE)
  stopifnot(all(c("gene_id", "cell_type") %in% colnames(module_markers)))
  module_lists <- split(module_markers$gene_id, module_markers$cell_type)
  module_correlations <- list()
  for (cell_type in names(module_lists)) {
    proportion_column <- safe_spotlight_name(cell_type)
    if (!proportion_column %in% colnames(visium_query[[]])) next
    features <- intersect(unique(module_lists[[cell_type]]),
                          rownames(visium_query[[query_assay]]))
    if (length(features) < 2L) next
    temporary_name <- paste0("MODULE_", make.names(cell_type), "_")
    visium_query <- AddModuleScore(
      visium_query, features = list(features), assay = query_assay,
      name = temporary_name, seed = 1
    )
    score_column <- paste0(temporary_name, "1")
    valid <- is.finite(visium_query[[proportion_column, drop = TRUE]]) &
      is.finite(visium_query[[score_column, drop = TRUE]])
    module_correlations[[cell_type]] <- data.frame(
      cell_type = cell_type,
      n_spots = sum(valid),
      spearman_rho = cor(
        visium_query[[proportion_column, drop = TRUE]][valid],
        visium_query[[score_column, drop = TRUE]][valid],
        method = "spearman"
      )
    )
  }
  write.csv(bind_rows(module_correlations),
            file.path(table_dir, "SPOTlight_module_score_correlations.csv"),
            row.names = FALSE)
}

# -----------------------------
# UMAP visualizations and original-study thresholds
# -----------------------------

query_reduction <- first_existing(
  c("umapharmony", "umap.harmony", "harmony.umap", "umap", "ref.umap"),
  Reductions(visium_query), "Visium UMAP reduction"
)
spotlight_columns <- grep("^SPOT_", colnames(visium_query[[]]), value = TRUE)
spotlight_columns <- setdiff(
  spotlight_columns,
  c("SPOT_top_type", "SPOT_top_prop", "SPOT_high_purity",
    "SPOT_entropy", "SPOT_entropy_normalized")
)
if (length(spotlight_columns)) {
  proportion_umap <- FeaturePlot(
    visium_query,
    features = spotlight_columns,
    reduction = query_reduction,
    order = TRUE,
    ncol = 4,
    min.cutoff = 0,
    max.cutoff = 1
  ) + plot_annotation(title = "SPOTlight cell-type proportions on Visium UMAP")
  ggsave(file.path(figure_dir, "SPOTlight_celltype_proportions_Visium_UMAP.png"),
         proportion_umap, width = 14,
         height = max(7, ceiling(length(spotlight_columns) / 4) * 3.2),
         dpi = 300)
}

plot_threshold_overlay <- function(object, feature, cutoff, label) {
  if (!feature %in% colnames(object[[]])) return(NULL)
  high_cells <- colnames(object)[
    is.finite(object[[feature, drop = TRUE]]) &
      object[[feature, drop = TRUE]] > cutoff
  ]
  DimPlot(
    object,
    reduction = query_reduction,
    cells.highlight = high_cells,
    cols = "grey82",
    cols.highlight = "#E41A1C",
    pt.size = 0.25,
    sizes.highlight = 0.45
  ) +
    NoLegend() +
    ggtitle(paste0(label, " (", feature, " > ", cutoff, ")"))
}

vascular_feature <- first_existing(
  c("SPOT_Vascular_tissue", "SPOT_Vascular"),
  colnames(visium_query[[]]), "vascular SPOTlight field", required = FALSE
)
sam_feature <- first_existing(
  c("SPOT_SAM_tissue", "SPOT_Shoot_apical_meristem", "SPOT_SAM"),
  colnames(visium_query[[]]), "SAM SPOTlight field", required = FALSE
)
threshold_plots <- Filter(Negate(is.null), list(
  if (!is.na(vascular_feature)) {
    plot_threshold_overlay(visium_query, vascular_feature,
                           vascular_display_cutoff, "Vascular-enriched spots")
  },
  if (!is.na(sam_feature)) {
    plot_threshold_overlay(visium_query, sam_feature,
                           sam_display_cutoff, "SAM-enriched spots")
  }
))
if (length(threshold_plots)) {
  threshold_figure <- wrap_plots(threshold_plots, ncol = 2)
  ggsave(file.path(figure_dir, "SPOTlight_original_threshold_UMAP_overlays.png"),
         threshold_figure, width = 11, height = 5.5, dpi = 300)
}

Idents(visium_query) <- original_visium_idents
stopifnot(identical(Idents(visium_query), original_visium_idents))
saveRDS(visium_query, output_rds, compress = FALSE)

writeLines(
  capture.output(sessionInfo()),
  file.path(session_dir, "12_scRNA_Visium_mapping_sessionInfo.txt")
)

message("scRNA-to-Visium mapping workflow completed.")
message("Output RDS: ", output_rds)
message("Results: ", figure_dir)
