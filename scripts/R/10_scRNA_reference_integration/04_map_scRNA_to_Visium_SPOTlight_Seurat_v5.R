# Map annotated scRNA-seq cell types onto Visium spots with Seurat and SPOTlight
#
# Main inputs:
#   data/processed/sc_merged_filter_SCT2_inte_SCINA.rds
#   data/processed/XGE202122_S5_subset_embleaf_harmony_join.rds
#
# Analysis scope:
#   Only the 6,392 embryonic-leaf spots assigned to SAM, P1_P2, P3, P4,
#   or P5 are used for both Seurat label transfer and SPOTlight deconvolution.
#   The `umap.harmony` coordinates stored in this subset are retained.
#
# Recommended for large objects:
#   Load both objects in RStudio as `sc_reference` and `visium_query` before
#   sourcing this script. In-memory objects take precedence over RDS files.
#
# Output:
#   data/processed/XGE202122_S5_subset_embleaf_celltype_mapped_SPOTlight_seurat_v5.rds

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
# FALSE intentionally recalculates transfer labels for the SAM-P5 subset rather
# than reusing labels previously calculated on the complete 14-sample object.
reuse_existing_hard_transfer <- FALSE
resume_spotlight_from_memory <- TRUE
# NULL uses the first available field in the documented preference order.
# Set this explicitly to "celltype_scina" for fully automated SCINA labels or
# "celltype_scina_histo" for the curated annotation used by the original
# 12-cell-type Figure B/C representation.
reference_celltype_column <- NULL
sections_to_run <- NULL       # NULL runs every section; use c("VR03_S2") to test
scatterpie_sections <- c("VR03_S1", "VR03_S2", "VR03_S3", "VR03_S4")
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
  gsub("^SPOT_", "", as_character_vector(x, "cell type"))
}

as_character_vector <- function(x, description = "value") {
  if (is.data.frame(x)) {
    if (ncol(x) != 1L) {
      stop(description, " must contain exactly one column.")
    }
    x <- x[[1L]]
  }
  if (methods::is(x, "Rle")) {
    # Factor-Rle objects do not consistently support as.character() across
    # S4Vectors versions, so expand their ordinary run values explicitly.
    x <- rep(
      as.character(S4Vectors::runValue(x)),
      S4Vectors::runLength(x)
    )
  }
  if (isS4(x)) {
    x <- tryCatch(
      as.vector(x),
      error = function(error) {
        stop(
          "Could not convert ", description, " from S4 class ",
          paste(class(x), collapse = "/"), " to an ordinary vector: ",
          conditionMessage(error)
        )
      }
    )
  }
  if (isS4(x) || is.list(x)) {
    stop(description, " is not an atomic vector after extraction.")
  }
  as.character(x)
}

update_legacy_visium_object <- function(object) {
  if (!inherits(object, "Seurat")) stop("Expected a Seurat object.")
  identities_before <- SeuratObject::Idents(object)
  identity_names <- names(identities_before)

  for (image_name in SeuratObject::Images(object)) {
    object@images[[image_name]] <- SeuratObject::UpdateSlots(
      object@images[[image_name]]
    )
  }
  object <- SeuratObject::UpdateSlots(object)
  object <- SeuratObject::UpdateSeuratObject(object)

  if (!is.null(identity_names) &&
      all(colnames(object) %in% identity_names)) {
    SeuratObject::Idents(object) <- identities_before[colnames(object)]
  } else {
    SeuratObject::Idents(object) <- identities_before
  }
  methods::validObject(object)
  object
}

get_all_tissue_coordinates <- function(object) {
  image_names <- SeuratObject::Images(object)
  if (!length(image_names)) stop("The Visium object contains no spatial images.")

  coordinate_list <- lapply(image_names, function(image_name) {
    coordinates <- SeuratObject::GetTissueCoordinates(
      object, image = image_name
    )
    if (!nrow(coordinates)) return(NULL)
    coordinates$spatial_image <- image_name
    coordinates$spot_barcode <- rownames(coordinates)
    coordinates
  })
  coordinates <- dplyr::bind_rows(coordinate_list)
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

build_marker_table <- function(sce, groups, genes_to_test,
                               maximum_per_type = 100L) {
  if (length(groups) != ncol(sce)) {
    stop("Marker-group vector length does not match the reference cells.")
  }
  if (!requireNamespace("presto", quietly = TRUE)) {
    stop(
      "The presto package is required for sparse AUC marker scoring. ",
      "Install presto before running SPOTlight reference preparation."
    )
  }
  message("Scoring reference markers by cell type with presto::wilcoxauc().")
  log_expression <- SummarizedExperiment::assay(
    sce, "logcounts", withDimnames = TRUE
  )[genes_to_test, , drop = FALSE]
  marker_statistics <- presto::wilcoxauc(
    log_expression,
    y = factor(groups),
    verbose = TRUE
  )
  message(
    "Marker scoring returned ", nrow(marker_statistics),
    " gene-by-cell-type statistics."
  )
  rm(log_expression)

  required_marker_columns <- c("feature", "group", "auc", "logFC")
  if (!all(required_marker_columns %in% colnames(marker_statistics))) {
    stop(
      "presto::wilcoxauc() did not return the expected columns: ",
      paste(required_marker_columns, collapse = ", ")
    )
  }
  marker_table <- marker_statistics |>
    dplyr::filter(
      feature %in% genes_to_test,
      is.finite(auc),
      is.finite(logFC),
      auc > 0.5,
      logFC > 0
    ) |>
    dplyr::arrange(group, dplyr::desc(auc), dplyr::desc(logFC)) |>
    dplyr::group_by(group) |>
    dplyr::slice_head(n = maximum_per_type) |>
    dplyr::ungroup() |>
    dplyr::transmute(
      gene = as.character(feature),
      cluster = as.character(group),
      weight = as.numeric(auc)
    )
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

build_serialization_safe_output <- function(mapped_object, base_rds,
                                            identities_to_restore) {
  if (!file.exists(base_rds)) {
    stop("Cannot construct the output because the base Visium RDS is missing: ",
         base_rds)
  }
  output_object <- readRDS(base_rds)
  if (!inherits(output_object, "Seurat")) {
    stop("The base Visium RDS does not contain a Seurat object.")
  }
  output_object <- update_legacy_visium_object(output_object)
  if (!setequal(colnames(output_object), colnames(mapped_object))) {
    stop("The base and mapped Visium objects do not contain the same spots.")
  }

  mapped_metadata <- mapped_object[[]]
  stale_output_columns <- grep(
    "^(predicted\\.|prediction\\.score\\.|mapping\\.score|SPOT_)",
    colnames(output_object[[]]), value = TRUE
  )
  if (length(stale_output_columns)) {
    for (field in stale_output_columns) output_object[[field]] <- NULL
  }
  output_columns <- grep(
    "^(predicted\\.|prediction\\.score\\.|mapping\\.score|SPOT_)",
    colnames(mapped_metadata),
    value = TRUE
  )
  if (!length(output_columns)) {
    stop("No transferred-label or SPOTlight metadata fields were found.")
  }
  output_metadata <- mapped_metadata[
    colnames(output_object), output_columns, drop = FALSE
  ]
  output_object <- SeuratObject::AddMetaData(
    output_object, metadata = output_metadata
  )

  identity_names <- names(identities_to_restore)
  if (!is.null(identity_names) &&
      all(colnames(output_object) %in% identity_names)) {
    SeuratObject::Idents(output_object) <- identities_to_restore[
      colnames(output_object)
    ]
  } else {
    SeuratObject::Idents(output_object) <- identities_to_restore
  }
  output_object
}

project_root <- find_project_root()
reference_rds <- file.path(
  project_root, "data", "processed",
  "sc_merged_filter_SCT2_inte_SCINA.rds"
)
visium_rds <- file.path(
  project_root, "data", "processed",
  "XGE202122_S5_subset_embleaf_harmony_join.rds"
)
output_rds <- file.path(
  project_root, "data", "processed",
  "XGE202122_S5_subset_embleaf_celltype_mapped_SPOTlight_seurat_v5.rds"
)
optional_module_marker_csv <- file.path(
  project_root, "data", "metadata", "scRNA_reference",
  "independent_celltype_module_markers.csv"
)
mitochondrial_gene_file <- file.path(
  project_root, "data", "reference", "maize_mitochondrial_genes.txt"
)
plastid_gene_file <- file.path(
  project_root, "data", "reference", "maize_plastid_genes.txt"
)
figure_dir <- file.path(
  project_root, "results", "figures", "10_scRNA_Visium_mapping"
)
table_dir <- file.path(
  project_root, "results", "tables", "10_scRNA_Visium_mapping"
)
proportion_checkpoint_rds <- file.path(
  table_dir, "SPOTlight_subset_embleaf_all_proportions_checkpoint.rds"
)
diagnostic_checkpoint_rds <- file.path(
  table_dir, "SPOTlight_subset_embleaf_section_diagnostics_checkpoint.rds"
)
session_dir <- file.path(project_root, "results", "sessionInfo")
section_pie_dir <- file.path(figure_dir, "section_scatterpies")

for (directory in c(dirname(output_rds), figure_dir, table_dir,
                    session_dir, section_pie_dir)) {
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
}

if (exists("sc_reference", envir = .GlobalEnv, inherits = FALSE)) {
  sc_reference <- get("sc_reference", envir = .GlobalEnv)
} else if (exists("sc_merged_filter_SCT2_inte_SCINA",
                  envir = .GlobalEnv, inherits = FALSE)) {
  sc_reference <- get(
    "sc_merged_filter_SCT2_inte_SCINA", envir = .GlobalEnv
  )
} else if (file.exists(reference_rds)) {
  sc_reference <- readRDS(reference_rds)
} else {
  stop(
    "The SCINA-annotated scRNA reference is required but was not found.\n",
    "First run:\n",
    "source(\"scripts/R/10_scRNA_reference_integration/",
    "03_scRNA_celltype_annotation_SCINA_Seurat_v5.R\")\n",
    "This creates: ", reference_rds, "\n",
    "Alternatively, load that annotated Seurat object as `sc_reference`."
  )
}

if (!file.exists(visium_rds)) {
  stop("The required SAM-P5 Visium subset was not found: ", visium_rds)
}

# The deposited subset is authoritative for both analysis membership and UMAP
# coordinates. An in-memory query may contain additional spots, so restrict it
# to the exact deposited subset before label transfer or deconvolution.
visium_scope <- readRDS(visium_rds)
stopifnot(inherits(visium_scope, "Seurat"))
scope_cells <- colnames(visium_scope)
scope_reduction <- "umap.harmony"
if (!scope_reduction %in% SeuratObject::Reductions(visium_scope)) {
  stop("The SAM-P5 subset lacks the required `umap.harmony` reduction.")
}
scope_embeddings <- SeuratObject::Embeddings(
  visium_scope, reduction = scope_reduction
)[, 1:2, drop = FALSE]

if (exists("visium_query", envir = .GlobalEnv, inherits = FALSE)) {
  visium_query <- get("visium_query", envir = .GlobalEnv)
  stopifnot(inherits(visium_query, "Seurat"))
  missing_scope_cells <- setdiff(scope_cells, colnames(visium_query))
  if (length(missing_scope_cells)) {
    stop(
      "The in-memory `visium_query` is missing ", length(missing_scope_cells),
      " spots required by the deposited SAM-P5 subset."
    )
  }
  visium_query <- subset(visium_query, cells = scope_cells)
} else {
  visium_query <- visium_scope
}

# Deposited legacy VisiumV1 images can lack the `misc` slot required by current
# SeuratObject. Upgrade all spatial images before MapQuery modifies the object.
visium_query <- update_legacy_visium_object(visium_query)

# Recreate the plotting reduction in the final query order so every UMAP plot
# uses the exact coordinates from XGE202122_S5_subset_embleaf_harmony_join.rds.
scope_embeddings <- scope_embeddings[colnames(visium_query), , drop = FALSE]
colnames(scope_embeddings) <- c("umapharmony_1", "umapharmony_2")
visium_query[[scope_reduction]] <- SeuratObject::CreateDimReducObject(
  embeddings = scope_embeddings,
  key = "umapharmony_",
  assay = SeuratObject::DefaultAssay(visium_query)
)
rm(visium_scope, scope_embeddings)

stopifnot(inherits(sc_reference, "Seurat"), inherits(visium_query, "Seurat"))
if (ncol(visium_query) != length(scope_cells) ||
    !setequal(colnames(visium_query), scope_cells)) {
  stop("The Visium query does not exactly match the deposited SAM-P5 subset.")
}

domain_column <- first_existing(
  c("domains", "domain", "structural_domain"), colnames(visium_query[[]]),
  "structural-domain metadata column"
)
allowed_domains <- c("SAM", "P1_P2", "P3", "P4", "P5")
query_domains <- as_character_vector(
  visium_query[[domain_column, drop = TRUE]], domain_column
)
unexpected_domains <- setdiff(
  unique(query_domains[!is.na(query_domains) & nzchar(query_domains)]),
  allowed_domains
)
if (length(unexpected_domains)) {
  stop(
    "The Visium query contains domains outside the SAM-P5 scope: ",
    paste(unexpected_domains, collapse = ", ")
  )
}
message(
  "Loaded annotated scRNA reference (", ncol(sc_reference),
  " cells) and SAM-P5 Visium query (", ncol(visium_query), " spots)."
)
original_visium_idents <- SeuratObject::Idents(visium_query)

# Remove deposited/previous mapping fields when the corresponding calculation
# is being redone. This prevents old full-query results from surviving in the
# subset output when a section is skipped or a field name changes.
if (run_hard_label_transfer && !reuse_existing_hard_transfer) {
  previous_transfer_columns <- grep(
    "^(predicted\\.|prediction\\.score\\.|mapping\\.score)",
    colnames(visium_query[[]]), value = TRUE
  )
  if (length(previous_transfer_columns)) {
    for (field in previous_transfer_columns) visium_query[[field]] <- NULL
  }
}
if (run_spotlight_deconvolution) {
  previous_spotlight_columns <- grep(
    "^SPOT_", colnames(visium_query[[]]), value = TRUE
  )
  if (length(previous_spotlight_columns)) {
    for (field in previous_spotlight_columns) visium_query[[field]] <- NULL
  }
}
reference_metadata <- sc_reference[[]]
visium_metadata <- visium_query[[]]

if (is.null(reference_celltype_column)) {
  celltype_column <- first_existing(
    c("celltype_scina", "celltype_scina_histo", "predicted_cell_type",
      "ident", "label"),
    colnames(reference_metadata), "scRNA cell-type annotation column"
  )
} else {
  if (!reference_celltype_column %in% colnames(reference_metadata)) {
    stop(
      "Requested reference cell-type column was not found: ",
      reference_celltype_column
    )
  }
  celltype_column <- reference_celltype_column
}
message("Using scRNA reference cell-type field: ", celltype_column)
reference_labels <- as_character_vector(
  reference_metadata[[celltype_column]], celltype_column
)
valid_reference_cells <- !is.na(reference_labels) & nzchar(reference_labels) &
  reference_labels != "Unknown"
sc_reference <- subset(sc_reference, cells = colnames(sc_reference)[valid_reference_cells])
sc_reference$celltype_mapping_reference <- factor(
  reference_labels[valid_reference_cells]
)

reference_count_assay <- first_existing(
  c("RNA", "SCT"), SeuratObject::Assays(sc_reference),
  "reference raw-count assay"
)
query_assay <- first_existing(
  c("SCT", "RNA", "Spatial"), SeuratObject::Assays(visium_query),
  "Visium expression assay"
)

# -----------------------------
# Hard label transfer with Seurat anchors and MapQuery
# -----------------------------

hard_transfer_available <- "predicted.celltype" %in%
  colnames(visium_query[[]])

if (run_hard_label_transfer && reuse_existing_hard_transfer &&
    hard_transfer_available) {
  message("Reusing the completed Seurat hard-label transfer in visium_query.")
} else if (run_hard_label_transfer) {
  message("Starting Seurat SCT anchor transfer using PCA projection.")
  if (!"SCT" %in% SeuratObject::Assays(sc_reference) ||
      !"SCT" %in% SeuratObject::Assays(visium_query)) {
    stop("SCT assays are required for the configured SCT anchor transfer.")
  }
  # Seurat's pcaproject mapping requires a PCA reduction with feature loadings.
  # Harmony is retained for visualization but is not used as the projection
  # basis for FindTransferAnchors/MapQuery.
  if (!"pca" %in% SeuratObject::Reductions(sc_reference)) {
    SeuratObject::DefaultAssay(sc_reference) <- "SCT"
    set.seed(1)
    sc_reference <- RunPCA(
      sc_reference, assay = "SCT", npcs = n_pcs,
      reduction.name = "pca", reduction.key = "PC_", verbose = FALSE
    )
  }
  reference_reduction <- "pca"
  available_dimensions <- min(
    n_pcs,
    ncol(SeuratObject::Embeddings(
      sc_reference, reduction = reference_reduction
    ))
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
    reduction = "pcaproject",
    reference.reduction = reference_reduction,
    dims = transfer_dimensions
  )

  message("Transfer anchors completed; projecting labels with MapQuery.")
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
  message("Starting SPOTlight reference preparation and deconvolution.")
  reference_cell_type_vector <- as_character_vector(
    sc_reference[[]][["celltype_mapping_reference"]],
    "scRNA reference cell types"
  )
  reference_counts <- SeuratObject::GetAssayData(
    sc_reference, assay = reference_count_assay, layer = "counts"
  )
  message(
    "Recovered raw reference counts from the ", reference_count_assay,
    " assay: ", nrow(reference_counts), " genes x ",
    ncol(reference_counts), " cells."
  )
  sce_reference <- SingleCellExperiment::SingleCellExperiment(
    assays = list(counts = reference_counts),
    colData = S4Vectors::DataFrame(
      cell_type = reference_cell_type_vector
    )
  )
  message("Created the SingleCellExperiment reference.")
  sce_reference <- scuttle::logNormCounts(sce_reference)
  message("Completed log-normalization for SPOTlight marker selection.")
  SingleCellExperiment::colLabels(sce_reference) <- factor(
    reference_cell_type_vector
  )
  message("Assigned reference cell-type labels.")

  mitochondrial_genes <- read_optional_gene_list(mitochondrial_gene_file)
  plastid_genes <- read_optional_gene_list(plastid_gene_file)
  eligible_genes <- remove_organelle_and_ribosomal_genes(
    rownames(sce_reference), mitochondrial_genes, plastid_genes
  )

  variance_model <- scran::modelGeneVar(
    sce_reference, subset.row = eligible_genes
  )
  message("Completed reference variance modeling.")
  requested_hvg <- min(n_hvg, nrow(variance_model))
  if (requireNamespace("scrapper", quietly = TRUE)) {
    hvg_indices <- scrapper::chooseHighlyVariableGenes(
      variance_model$bio,
      top = requested_hvg,
      bound = 0
    )
    hvg <- rownames(variance_model)[hvg_indices]
  } else {
    hvg <- scran::getTopHVGs(variance_model, n = requested_hvg)
  }
  message("Selected ", length(hvg), " highly variable genes.")
  marker_table <- build_marker_table(
    sce_reference,
    groups = reference_cell_type_vector,
    genes_to_test = eligible_genes,
    maximum_per_type = markers_per_type
  )
  message("Completed SPOTlight marker and HVG selection.")
  write.csv(marker_table,
            file.path(table_dir, "SPOTlight_scRNA_marker_table.csv"),
            row.names = FALSE)
  write.csv(data.frame(gene = hvg),
            file.path(table_dir, "SPOTlight_scRNA_HVGs.csv"),
            row.names = FALSE)

  sce_reference_small <- downsample_reference(
    sce_reference,
    groups = as_character_vector(
      SummarizedExperiment::colData(sce_reference)$cell_type,
      "scRNA reference cell types"
    ),
    maximum_per_type = reference_cells_per_type
  )
  reference_cell_types <- levels(droplevels(
    factor(as_character_vector(
      SummarizedExperiment::colData(sce_reference_small)$cell_type,
      "downsampled scRNA reference cell types"
    ))
  ))

  spatial_assay <- first_existing(
    c("RNA", "Spatial"), SeuratObject::Assays(visium_query),
    "raw Visium count assay"
  )
  spatial_count_layers <- grep(
    "^counts($|\\.)",
    SeuratObject::Layers(visium_query[[spatial_assay]]),
    value = TRUE
  )
  if (!length(spatial_count_layers)) {
    stop("No raw count layer was found in the Visium ", spatial_assay, " assay.")
  }
  if (length(spatial_count_layers) > 1L) {
    message(
      "Joining ", length(spatial_count_layers), " Visium ", spatial_assay,
      " count layers for SPOTlight."
    )
    visium_query[[spatial_assay]] <- SeuratObject::JoinLayers(
      visium_query[[spatial_assay]],
      # JoinLayers expects the layer-family search term, not the expanded list
      # of counts.sample layer names, in SeuratObject 5.4.
      layers = "counts",
      new = "counts"
    )
    joined_count_layers <- grep(
      "^counts($|\\.)",
      SeuratObject::Layers(visium_query[[spatial_assay]]),
      value = TRUE
    )
    if (!identical(joined_count_layers, "counts")) {
      stop(
        "Visium count-layer joining did not produce one `counts` layer. ",
        "Remaining layers: ", paste(joined_count_layers, collapse = ", ")
      )
    }
  }
  spatial_counts <- SeuratObject::GetAssayData(
    visium_query, assay = spatial_assay, layer = "counts"
  )
  message(
    "Recovered joined Visium counts: ", nrow(spatial_counts),
    " genes x ", ncol(spatial_counts), " spots."
  )
  all_coordinates <- get_all_tissue_coordinates(visium_query)
  xy_columns <- select_xy_columns(all_coordinates)

  section_column <- first_existing(
    c("section_id", "domain_section", "section"),
    colnames(visium_query[[]]), "section metadata column"
  )
  section_labels <- as_character_vector(
    visium_query[[section_column, drop = TRUE]], section_column
  )
  names(section_labels) <- colnames(visium_query)
  available_sections <- sort(unique(section_labels[!is.na(section_labels) &
                                                      nzchar(section_labels)]))
  target_sections <- if (is.null(sections_to_run)) {
    available_sections
  } else {
    intersect(sections_to_run, available_sections)
  }
  if (!length(target_sections)) stop("No requested Visium sections were found.")

  previous_proportions <- if (
    resume_spotlight_from_memory &&
      exists("all_proportions", envir = .GlobalEnv, inherits = FALSE)
  ) {
    get("all_proportions", envir = .GlobalEnv)
  } else if (file.exists(proportion_checkpoint_rds)) {
    readRDS(proportion_checkpoint_rds)
  } else {
    NULL
  }
  valid_previous_proportions <- is.matrix(previous_proportions) &&
    identical(rownames(previous_proportions), colnames(visium_query)) &&
    identical(colnames(previous_proportions), reference_cell_types)

  if (valid_previous_proportions) {
    all_proportions <- previous_proportions
    message("Recovered an existing SPOTlight proportion checkpoint.")
  } else {
    all_proportions <- matrix(
      NA_real_, nrow = ncol(visium_query), ncol = length(reference_cell_types),
      dimnames = list(colnames(visium_query), reference_cell_types)
    )
  }

  previous_diagnostics <- if (
    resume_spotlight_from_memory &&
      exists("section_diagnostics", envir = .GlobalEnv, inherits = FALSE)
  ) {
    get("section_diagnostics", envir = .GlobalEnv)
  } else if (file.exists(diagnostic_checkpoint_rds)) {
    readRDS(diagnostic_checkpoint_rds)
  } else {
    NULL
  }
  if (is.list(previous_diagnostics) &&
      all(target_sections %in% names(previous_diagnostics))) {
    section_diagnostics <- previous_diagnostics[target_sections]
  } else {
    section_diagnostics <- vector("list", length(target_sections))
    names(section_diagnostics) <- target_sections
  }

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
    existing_section <- all_proportions[section_spots, , drop = FALSE]
    if (all(rowSums(existing_section, na.rm = TRUE) > 0)) {
      message("Reusing completed SPOTlight proportions for: ", section_id)
      next
    }

    common_genes <- intersect(
      rownames(sce_reference_small), rownames(spatial_counts)
    )
    section_markers <- marker_table |>
      dplyr::filter(gene %in% common_genes)
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
        as_character_vector(
          SummarizedExperiment::colData(section_reference)$cell_type,
          "section reference cell types"
        )
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

    # Persist section-level progress so an interrupted R session can resume
    # without repeating completed NMF/deconvolution calculations.
    saveRDS(all_proportions, proportion_checkpoint_rds, compress = FALSE)
    saveRDS(section_diagnostics, diagnostic_checkpoint_rds, compress = FALSE)
    message("Saved SPOTlight checkpoint after: ", section_id)

  }

  message("All requested SPOTlight section calculations are complete.")
  completed_spots <- rowSums(all_proportions, na.rm = TRUE) > 0
  proportion_matrix <- all_proportions[completed_spots, , drop = FALSE]
  if (!nrow(proportion_matrix)) {
    stop("SPOTlight did not return a positive proportion vector for any spot.")
  }
  proportion_matrix <- normalize_rows(proportion_matrix)
  all_proportions[completed_spots, ] <- proportion_matrix
  saveRDS(all_proportions, proportion_checkpoint_rds, compress = FALSE)
  saveRDS(section_diagnostics, diagnostic_checkpoint_rds, compress = FALSE)
  message(
    "Saved the complete SPOTlight proportion checkpoint for ",
    sum(completed_spots), " spots."
  )
  proportion_metadata <- all_proportions
  colnames(proportion_metadata) <- safe_spotlight_name(
    colnames(proportion_metadata)
  )
  message("Attaching SPOTlight proportions to the Visium metadata.")
  visium_query <- SeuratObject::AddMetaData(
    visium_query, proportion_metadata
  )
  message("Attached SPOTlight proportion fields.")

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
  visium_query <- SeuratObject::AddMetaData(visium_query, summary_metadata)
  message("Attached SPOTlight summary fields.")

  write.csv(
    cbind(spot = rownames(proportion_matrix), as.data.frame(proportion_matrix)),
    file.path(table_dir, "SPOTlight_subset_embleaf_spot_celltype_proportions.csv"),
    row.names = FALSE
  )
  write.csv(
    dplyr::bind_rows(section_diagnostics),
    file.path(table_dir, "SPOTlight_subset_embleaf_section_diagnostics.csv"),
    row.names = FALSE
  )

  # MapQuery can add transient projected reductions whose embedded UMAP model
  # is not serializable in some Seurat/uwot combinations. Transfer the durable
  # mapped metadata onto a clean copy of the original, already serializable
  # Visium object. Its assays, images, Harmony reduction, and identities remain
  # unchanged; only predicted-label and SPOTlight fields are added.
  message("Building a serialization-safe mapped Visium object.")
  visium_query <- build_serialization_safe_output(
    mapped_object = visium_query,
    base_rds = visium_rds,
    identities_to_restore = original_visium_idents
  )
  message("Saving the SPOTlight-mapped Visium checkpoint.")
  saveRDS(visium_query, output_rds, compress = FALSE)
  message("Saved the SPOTlight-mapped Visium checkpoint: ", output_rds)

  # Generate only the four serial VR03 scatter-pie panels used for Figure B.
  sections_for_pies <- intersect(scatterpie_sections, target_sections)
  for (section_id in sections_for_pies) {
    section_spots <- names(section_labels)[section_labels == section_id]
    section_spots <- intersect(
      section_spots,
      intersect(rownames(all_coordinates), rownames(all_proportions))
    )
    coordinate_data <- all_coordinates[section_spots, xy_columns, drop = FALSE]
    coordinate_matrix <- cbind(
      x = suppressWarnings(as.numeric(coordinate_data[[xy_columns[1L]]])),
      y = suppressWarnings(as.numeric(coordinate_data[[xy_columns[2L]]]))
    )
    rownames(coordinate_matrix) <- section_spots
    pie_proportions <- all_proportions[section_spots, , drop = FALSE]
    valid_plot_spots <- rowSums(is.finite(coordinate_matrix)) == 2L &
      rowSums(pie_proportions, na.rm = TRUE) > 0
    coordinate_matrix <- coordinate_matrix[valid_plot_spots, , drop = FALSE]
    pie_proportions <- pie_proportions[valid_plot_spots, , drop = FALSE]

    if (!nrow(pie_proportions)) {
      warning("Skipping scatter-pie plot for ", section_id,
              ": no finite aligned coordinates.")
      next
    }
    tryCatch({
      scatterpie_plot <- SPOTlight::plotSpatialScatterpie(
        x = coordinate_matrix,
        y = pie_proportions,
        cell_types = colnames(pie_proportions),
        img = FALSE,
        scatterpie_alpha = 0.85,
        pie_scale = 0.45
      ) + ggtitle(paste("SPOTlight cell-type proportions:", section_id))
      ggsave(
        file.path(section_pie_dir, paste0(make.names(section_id), ".png")),
        scatterpie_plot, width = 8, height = 7, dpi = 300
      )
    }, error = function(error) {
      warning("Scatter-pie plotting failed for ", section_id, ": ",
              conditionMessage(error))
    })
  }
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
            file.path(
              table_dir,
              "SPOTlight_subset_embleaf_vs_Seurat_label_agreement.csv"
            ),
            row.names = FALSE)
  write.csv(
    data.frame(
      n_evaluable_spots = nrow(agreement_table),
      agreement_rate = mean(agreement_table$agreement)
    ),
    file.path(
      table_dir,
      "SPOTlight_subset_embleaf_vs_Seurat_agreement_summary.csv"
    ),
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
  write.csv(dplyr::bind_rows(module_correlations),
            file.path(
              table_dir,
              "SPOTlight_subset_embleaf_module_score_correlations.csv"
            ),
            row.names = FALSE)
}

# -----------------------------
# UMAP visualizations and original-study thresholds
# -----------------------------

query_reduction <- "umap.harmony"
if (!query_reduction %in% SeuratObject::Reductions(visium_query)) {
  stop("The mapped SAM-P5 query lacks the required `umap.harmony` reduction.")
}
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

SeuratObject::Idents(visium_query) <- original_visium_idents
stopifnot(identical(
  SeuratObject::Idents(visium_query), original_visium_idents
))
saveRDS(visium_query, output_rds, compress = FALSE)

writeLines(
  capture.output(sessionInfo()),
  file.path(session_dir, "12_scRNA_Visium_mapping_sessionInfo.txt")
)

message("scRNA-to-Visium mapping workflow completed.")
message("Output RDS: ", output_rds)
message("Results: ", figure_dir)
