#!/usr/bin/env Rscript

# Merge 14 individual maize shoot Seurat v5 objects and perform
# SCTransform normalization followed by Harmony integration, graph construction,
# and unsupervised clustering.
# After merging, cell names are standardized as BARCODE-1_1_SAMPLE_NUMBER
# (for example, UL01_AAACAAGTATCTCCCA-1 becomes
# AAACAAGTATCTCCCA-1_1_1). The original merged name is retained in
# meta.data$old_colname.
#
# Expected file structure (run this script from the repository root):
#
# maize_shoot_data_process_v2/
# ├── data/
# │   └── processed/
# │       ├── UL01_seurat_v5.rds
# │       ├── UL02_seurat_v5.rds
# │       ├── UL04_seurat_v5.rds
# │       ├── VR01_seurat_v5.rds
# │       ├── VR02_seurat_v5.rds
# │       ├── VR03_seurat_v5.rds
# │       ├── VR04_seurat_v5.rds
# │       ├── DQ01_seurat_v5.rds
# │       ├── DQ02_seurat_v5.rds
# │       ├── DQ03_seurat_v5.rds
# │       ├── DQ04_seurat_v5.rds
# │       ├── DQ06_seurat_v5.rds
# │       ├── DQ07_seurat_v5.rds
# │       ├── DQ08_seurat_v5.rds
# │       └── maize_shoot_14samples_SCT_harmony_seurat_v5.rds
# ├── data/metadata/metadata.csv
# └── scripts/
#     └── R/
#         └── 04_merge_integration/
#             └── 01_merge_14_samples_SCT_Harmony_Seurat_v5.R
#
# UL03 and DQ05 are not included because individual RDS files were not
# provided for these capture areas.
#
# Main outputs:
#   data/processed/maize_shoot_14samples_SCT_harmony_seurat_v5.rds
#   results/figures/PCA_elbow_plot.png
#   results/figures/PCA_Harmony_before_after_UMAP.png
#   results/figures/Harmony_UMAP_by_domain.png
#   results/figures/Harmony_UMAP_domains_by_sample.png
#   results/tables/harmony_clusters_recomputed_vs_published.csv
#   results/logs/Harmony_recomputed_clustering.txt

suppressPackageStartupMessages({
    library(Seurat)
    library(harmony)
    library(Matrix)
    library(patchwork)
})

# Test 50 PCs in the elbow plot and retain the first 30 PCs, matching the
# Bio-protocol and original analysis. The automatic elbow is also recorded as
# a diagnostic so that users can reassess this choice for other datasets.
max_pcs_to_test <- 50L
n_pcs_use <- 30L
minimum_total_reads_per_gene <- 100
random_seed <- 1234L
neighbor_k <- 20L
clustering_resolution <- 2
clustering_algorithm <- 1L

sample_ids <- c(
    "UL01", "UL02", "UL04",
    "VR01", "VR02", "VR03", "VR04",
    "DQ01", "DQ02", "DQ03", "DQ04", "DQ06", "DQ07", "DQ08"
)

allowed_domains <- c(
    "SAM", "P1_P2", "P3", "P4", "P5", "coleoptile", "co_v"
)
sample_domain_levels <- unlist(
    lapply(sample_ids, function(sample_id) {
        paste(sample_id, allowed_domains, sep = "_")
    }),
    use.names = FALSE
)
domain_colours <- setNames(
    c(
        "#1B9E77", "#D95F02", "#7570B3", "#E7298A",
        "#66A61E", "#E6AB02", "#A6761D"
    ),
    allowed_domains
)

rds_files <- file.path(
    "data",
    "processed",
    paste0(sample_ids, "_seurat_v5.rds")
)
domain_metadata <- read.csv(
    file.path("data", "metadata", "metadata.csv"),
    check.names = FALSE,
    stringsAsFactors = FALSE
)

stopifnot(all(c("Barcode", "domains") %in% colnames(domain_metadata)))
domain_metadata$Barcode <- trimws(domain_metadata$Barcode)
domain_metadata$domains <- trimws(domain_metadata$domains)
domain_metadata <- domain_metadata[
    domain_metadata$domains %in% allowed_domains,
    ,
    drop = FALSE
]
stopifnot(
    nrow(domain_metadata) > 0L,
    !anyNA(domain_metadata$Barcode),
    !anyDuplicated(domain_metadata$Barcode)
)
mitochondrial_genes <- scan(
    file.path("data", "reference", "maize_mitochondrial_genes.txt"),
    what = character(), quiet = TRUE
)
plastid_genes <- scan(
    file.path("data", "reference", "maize_plastid_genes.txt"),
    what = character(), quiet = TRUE
)

# Load the 14 individual objects and retain sample identity in metadata.
seurat_list <- lapply(seq_along(sample_ids), function(i) {
    object <- readRDS(rds_files[i])
    object <- UpdateSeuratObject(object)

    # Standardize older objects that used the default Spatial assay name.
    # SummarizedExperiment also exports Assays(), and it can mask the Seurat
    # method in an interactive session. Namespace this call so that it always
    # returns the character vector of Seurat assay names.
    assay_names <- SeuratObject::Assays(object)
    if (!"RNA" %in% assay_names && "Spatial" %in% assay_names) {
        object <- RenameAssays(object, Spatial = "RNA")
    }

    object$sample_id <- sample_ids[i]
    object$orig.ident <- sample_ids[i]
    DefaultAssay(object) <- "RNA"

    mt_use <- intersect(mitochondrial_genes, rownames(object[["RNA"]]))
    pt_use <- intersect(plastid_genes, rownames(object[["RNA"]]))
    stopifnot(length(mt_use) > 0L, length(pt_use) > 0L)
    object[["percent.mito"]] <- PercentageFeatureSet(
        object, assay = "RNA", features = mt_use
    )
    object[["percent.pltd"]] <- PercentageFeatureSet(
        object, assay = "RNA", features = pt_use
    )
    object
})
names(seurat_list) <- sample_ids

# In Seurat v5, merge() retains the individual datasets as separate RNA layers.
combined <- merge(
    x = seurat_list[[1]],
    y = seurat_list[-1],
    add.cell.ids = sample_ids,
    project = "maize_shoot_14samples",
    merge.data = FALSE,
    merge.dr = FALSE
)

# Standardize cell names for the combined 14-sample object. The sample order in
# sample_ids defines the final suffix: UL01 = _1_1, UL02 = _1_2, ...,
# DQ08 = _1_14. Preserve the original SAMPLE_BARCODE name in metadata.
old_cell_names <- colnames(combined)
cell_sample <- sub("_.*$", "", old_cell_names)
cell_barcode <- sub("^[^_]+_", "", old_cell_names)
cell_sample_number <- match(cell_sample, sample_ids)

if (anyNA(cell_sample_number)) {
    unknown_samples <- unique(cell_sample[is.na(cell_sample_number)])
    stop(
        "Unknown sample prefix in merged cell names: ",
        paste(unknown_samples, collapse = ", ")
    )
}

new_cell_names <- paste0(
    cell_barcode,
    "_1_",
    cell_sample_number
)

if (anyDuplicated(new_cell_names)) {
    stop("The standardized cell names are not unique.")
}

combined$old_colname <- old_cell_names
combined <- RenameCells(combined, new.names = new_cell_names)

stopifnot(
    identical(colnames(combined), new_cell_names),
    identical(
        unname(as.character(combined$old_colname)),
        unname(old_cell_names)
    ),
    identical(rownames(combined[[]]), new_cell_names)
)

# Retain only the non-overlapping shoot spots listed in metadata.csv, then align
# and add every metadata column by Barcode.
missing_metadata_spots <- setdiff(domain_metadata$Barcode, colnames(combined))
if (length(missing_metadata_spots) > 0L) {
    stop(
        length(missing_metadata_spots),
        " metadata barcodes are absent from the merged Seurat object."
    )
}

n_spots_before_domain_filter <- ncol(combined)
spots_to_keep <- colnames(combined)[
    colnames(combined) %in% domain_metadata$Barcode
]
combined <- subset(combined, cells = spots_to_keep)

metadata_aligned <- domain_metadata[
    match(colnames(combined), domain_metadata$Barcode),
    ,
    drop = FALSE
]
stopifnot(identical(metadata_aligned$Barcode, colnames(combined)))

# Preserve current pipeline-derived metadata under its canonical name. Historical
# CSV values for these fields are added with a _metadata_csv suffix.
protected_current_columns <- c(
    "orig.ident", "nCount_RNA", "nFeature_RNA",
    "nCount_SCT", "nFeature_SCT"
)
metadata_target_columns <- ifelse(
    colnames(metadata_aligned) %in% protected_current_columns,
    paste0(colnames(metadata_aligned), "_metadata_csv"),
    colnames(metadata_aligned)
)

for (i in seq_along(metadata_target_columns)) {
    combined[[metadata_target_columns[i]]] <- metadata_aligned[[i]]
}
combined$domains <- factor(
    as.character(combined$domains),
    levels = allowed_domains
)

fixed_level_factor_columns <- c("sample_id", "sample_domain")
missing_fixed_level_factor_columns <- setdiff(
    fixed_level_factor_columns,
    colnames(combined[[]])
)
if (length(missing_fixed_level_factor_columns) > 0L) {
    stop(
        "Metadata columns requested as fixed-level factors are missing: ",
        paste(missing_fixed_level_factor_columns, collapse = ", ")
    )
}
combined$sample_id <- factor(
    as.character(combined$sample_id),
    levels = sample_ids
)
combined$sample_domain <- factor(
    as.character(combined$sample_domain),
    levels = sample_domain_levels
)

factor_metadata_columns <- c(
    "orig.ident_metadata_csv", "section_id", "section", "cca_clusters",
    "seurat_clusters", "harmony_clusters", "ms_ve", "domain_section"
)
missing_factor_columns <- setdiff(
    factor_metadata_columns,
    colnames(combined[[]])
)
if (length(missing_factor_columns) > 0L) {
    stop(
        "Metadata columns requested as factors are missing: ",
        paste(missing_factor_columns, collapse = ", ")
    )
}
for (metadata_column in factor_metadata_columns) {
    combined[[metadata_column]] <- factor(
        combined[[metadata_column, drop = TRUE]]
    )
}

stopifnot(
    ncol(combined) == nrow(domain_metadata),
    !anyNA(combined$domains),
    all(as.character(combined$domains) %in% allowed_domains),
    is.factor(combined$sample_id),
    is.factor(combined$sample_domain),
    !anyNA(combined$sample_id),
    !anyNA(combined$sample_domain),
    identical(levels(combined$sample_id), sample_ids),
    identical(levels(combined$sample_domain), sample_domain_levels),
    all(metadata_target_columns %in% colnames(combined[[]])),
    all(vapply(
        factor_metadata_columns,
        function(x) is.factor(combined[[x, drop = TRUE]]),
        logical(1)
    ))
)
message(
    "Domain filtering retained ", ncol(combined), " of ",
    n_spots_before_domain_filter, " spots; removed ",
    n_spots_before_domain_filter - ncol(combined),
    " overlapping or non-shoot spots. Added all ",
    ncol(domain_metadata), " metadata.csv columns."
)

# Remove the temporary list after merging to avoid retaining duplicate objects.
rm(seurat_list)
gc()

DefaultAssay(combined) <- "RNA"

# The Bio-protocol removes genes with fewer than 100 reads across the combined
# dataset. Sum the split Seurat v5 count layers without joining or duplicating
# the complete sparse count matrix.
all_features <- rownames(combined[["RNA"]])
gene_total_reads <- setNames(numeric(length(all_features)), all_features)
count_layers <- Layers(combined[["RNA"]], search = "^counts")
for (layer_name in count_layers) {
    counts <- LayerData(combined, assay = "RNA", layer = layer_name)
    layer_features <- rownames(counts)
    gene_total_reads[layer_features] <- gene_total_reads[layer_features] +
        Matrix::rowSums(counts)
}
features_to_keep <- names(gene_total_reads)[
    gene_total_reads >= minimum_total_reads_per_gene
]
stopifnot(length(features_to_keep) > 0L)
combined <- subset(combined, features = features_to_keep)
message(
    "Retained ", length(features_to_keep), " genes with at least ",
    minimum_total_reads_per_gene, " total reads."
)

# Normalize the split sample layers using SCTransform v2.
combined <- SCTransform(
    object = combined,
    assay = "RNA",
    new.assay.name = "SCT",
    vst.flavor = "v2",
    variable.features.n = 3000,
    verbose = TRUE
)

# Calculate a broad set of PCs so that the appropriate number can be assessed.
combined <- RunPCA(
    object = combined,
    assay = "SCT",
    npcs = max_pcs_to_test,
    verbose = TRUE
)

# Summarize the variance explained by each PC.
pca_sd <- Stdev(combined, reduction = "pca")
pca_variance <- pca_sd^2
pca_percent <- 100 * pca_variance / sum(pca_variance)
pca_table <- data.frame(
    PC = seq_along(pca_sd),
    standard_deviation = pca_sd,
    variance_percent = pca_percent,
    cumulative_variance_percent = cumsum(pca_percent)
)

# Detect the elbow as the point with the greatest perpendicular distance from
# the straight line joining the first and last points of the normalized scree
# curve. This provides a reproducible candidate that is also shown on the plot.
x_scaled <- (pca_table$PC - min(pca_table$PC)) /
    (max(pca_table$PC) - min(pca_table$PC))
y_scaled <- (pca_sd - min(pca_sd)) / (max(pca_sd) - min(pca_sd))
elbow_distance <- abs(x_scaled + y_scaled - 1) / sqrt(2)
automatic_elbow_pc <- which.max(elbow_distance)

stopifnot(n_pcs_use >= 2L, n_pcs_use <= max_pcs_to_test)
selection_method <- "Bio-protocol setting"

dir.create(file.path("results", "tables"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path("results", "figures"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path("results", "logs"), recursive = TRUE, showWarnings = FALSE)

write.csv(
    pca_table,
    file = file.path("results", "tables", "PCA_variance_explained.csv"),
    row.names = FALSE
)

elbow_plot <- ElbowPlot(combined, ndims = max_pcs_to_test) +
    ggplot2::geom_vline(
        xintercept = automatic_elbow_pc,
        colour = "#2C7FB8",
        linetype = "dotted"
    ) +
    ggplot2::geom_vline(
        xintercept = n_pcs_use,
        colour = "red",
        linetype = "dashed"
    ) +
    ggplot2::labs(
        title = "PCA elbow plot",
        subtitle = paste0(
            "Retained PCs: ", n_pcs_use,
            " (Bio-protocol); automatic elbow candidate: ",
            automatic_elbow_pc
        )
    )

ggplot2::ggsave(
    filename = file.path("results", "figures", "PCA_elbow_plot.png"),
    plot = elbow_plot,
    width = 7,
    height = 5,
    dpi = 300,
    bg = "white"
)

writeLines(
    c(
        paste0("PCs tested: 1-", max_pcs_to_test),
        paste0("Automatic elbow candidate: ", automatic_elbow_pc),
        paste0("PCs used for Harmony: ", n_pcs_use),
        paste0("Selection method: ", selection_method)
    ),
    con = file.path("results", "logs", "PCA_selection.txt")
)

message("Automatic elbow candidate: PC ", automatic_elbow_pc)
message("Using PCs 1-", n_pcs_use, " for Harmony (", selection_method, ").")

# Retain only the selected PCs in the final object.
combined <- RunPCA(
    object = combined,
    assay = "SCT",
    npcs = n_pcs_use,
    verbose = TRUE
)

# UMAP before integration, using the uncorrected PCA space.
set.seed(random_seed)
combined <- RunUMAP(
    object = combined,
    reduction = "pca",
    dims = seq_len(n_pcs_use),
    n.neighbors = 30L,
    min.dist = 0.3,
    metric = "cosine",
    reduction.name = "umap_pca",
    reduction.key = "UMAPPCA_",
    seed.use = random_seed,
    verbose = FALSE
)

# Seurat v5 streamlined one-line Harmony integration.
combined <- IntegrateLayers(
    object = combined,
    method = HarmonyIntegration,
    orig.reduction = "pca",
    new.reduction = "harmony",
    assay = "SCT",
    dims = seq_len(n_pcs_use),
    verbose = FALSE
)

# Construct the shared-nearest-neighbor graph and calculate clusters from the
# newly generated Harmony reduction. Store these clusters under a new metadata
# name so the published demonstration labels imported from metadata.csv remain
# unchanged. Resolution 2 matches the original exploratory workflow, but users
# should evaluate an appropriate resolution for their own dataset.
published_harmony_clusters_before_reclustering <- as.character(
    combined$harmony_clusters
)
published_seurat_clusters_before_reclustering <- combined$seurat_clusters
combined <- FindNeighbors(
    object = combined,
    reduction = "harmony",
    dims = seq_len(n_pcs_use),
    k.param = neighbor_k,
    graph.name = c("harmony_nn", "harmony_snn"),
    verbose = TRUE
)
set.seed(random_seed)
combined <- FindClusters(
    object = combined,
    graph.name = "harmony_snn",
    resolution = clustering_resolution,
    algorithm = clustering_algorithm,
    random.seed = random_seed,
    cluster.name = "harmony_clusters_recomputed",
    verbose = TRUE
)
recomputed_cluster_values <- as.character(
    combined$harmony_clusters_recomputed
)
recomputed_cluster_levels <- as.character(sort(unique(as.integer(
    recomputed_cluster_values
))))
combined$harmony_clusters_recomputed <- factor(
    recomputed_cluster_values,
    levels = recomputed_cluster_levels
)
# FindClusters() also refreshes the generic seurat_clusters field. Restore the
# published value imported from metadata.csv so historical reference metadata
# remain unchanged; the new clustering is retained under its explicit name.
combined$seurat_clusters <- published_seurat_clusters_before_reclustering
stopifnot(
    identical(
        as.character(combined$harmony_clusters),
        published_harmony_clusters_before_reclustering
    ),
    identical(
        as.character(combined$seurat_clusters),
        as.character(published_seurat_clusters_before_reclustering)
    ),
    !anyNA(combined$harmony_clusters_recomputed)
)

# Record correspondence without assuming that cluster numbers are transferable.
# The recomputed labels are the appropriate starting point for a new analysis;
# the published labels are retained only to reproduce the manuscript example.
cluster_comparison <- as.data.frame(
    table(
        harmony_clusters_recomputed = combined$harmony_clusters_recomputed,
        harmony_clusters_published = combined$harmony_clusters
    ),
    stringsAsFactors = FALSE
)
write.csv(
    cluster_comparison,
    file = file.path(
        "results",
        "tables",
        "harmony_clusters_recomputed_vs_published.csv"
    ),
    row.names = FALSE
)
writeLines(
    c(
        paste0("Reduction: harmony"),
        paste0("Dimensions: 1-", n_pcs_use),
        paste0("k.param: ", neighbor_k),
        paste0("Graph: harmony_snn"),
        paste0("Clustering resolution: ", clustering_resolution),
        paste0("Clustering algorithm: ", clustering_algorithm),
        paste0("Random seed: ", random_seed),
        paste0(
            "Recomputed cluster count: ",
            nlevels(combined$harmony_clusters_recomputed)
        ),
        paste0(
            "Published cluster count: ",
            nlevels(combined$harmony_clusters)
        ),
        "harmony_clusters_recomputed was calculated from the current Harmony graph.",
        "harmony_clusters was imported from metadata.csv for manuscript reproduction."
    ),
    con = file.path(
        "results",
        "logs",
        "Harmony_recomputed_clustering.txt"
    )
)

# UMAP after SCTransform + Harmony integration.
set.seed(random_seed)
combined <- RunUMAP(
    object = combined,
    reduction = "harmony",
    dims = seq_len(n_pcs_use),
    n.neighbors = 30L,
    min.dist = 0.3,
    metric = "cosine",
    reduction.name = "umap_harmony",
    reduction.key = "UMAPHARMONY_",
    seed.use = random_seed,
    verbose = FALSE
)

before_integration <- DimPlot(
    combined,
    reduction = "umap_pca",
    group.by = "sample_id",
    pt.size = 0.2
) +
    ggplot2::ggtitle("Before integration") +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"))

after_integration <- DimPlot(
    combined,
    reduction = "umap_harmony",
    group.by = "sample_id",
    pt.size = 0.2
) +
    ggplot2::ggtitle("After SCTransform + Harmony") +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"))

integration_plot <- before_integration | after_integration
integration_plot <- integration_plot + plot_annotation(tag_levels = "A")
ggplot2::ggsave(
    filename = file.path(
        "results", "figures", "PCA_Harmony_before_after_UMAP.png"
    ),
    plot = integration_plot,
    width = 13,
    height = 6,
    dpi = 300,
    bg = "white"
)

# Confirm that Harmony preserves interpretable anatomical structure. A global
# domain plot shows the shared embedding, while the split view reveals whether
# the same domains occupy comparable regions in each capture area.
harmony_by_domain <- DimPlot(
    combined,
    reduction = "umap_harmony",
    group.by = "domains",
    cols = domain_colours,
    order = allowed_domains,
    label = TRUE,
    repel = TRUE,
    pt.size = 0.2,
    shuffle = TRUE,
    seed = random_seed
) +
    ggplot2::ggtitle("Harmony UMAP by anatomical domain") +
    ggplot2::theme(
        plot.title = ggplot2::element_text(hjust = 0.5, face = "bold")
    )
ggplot2::ggsave(
    filename = file.path("results", "figures", "Harmony_UMAP_by_domain.png"),
    plot = harmony_by_domain,
    width = 8,
    height = 6,
    dpi = 300,
    bg = "white"
)

harmony_domains_by_sample <- DimPlot(
    combined,
    reduction = "umap_harmony",
    group.by = "domains",
    split.by = "sample_id",
    cols = domain_colours,
    order = allowed_domains,
    ncol = 4,
    pt.size = 0.15,
    shuffle = TRUE,
    seed = random_seed
)
harmony_domains_by_sample <- (
    harmony_domains_by_sample +
        plot_annotation(title = "Harmony UMAP domains in each sample")
) &
    ggplot2::theme(
        plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
        legend.position = "bottom"
    )
ggplot2::ggsave(
    filename = file.path(
        "results", "figures", "Harmony_UMAP_domains_by_sample.png"
    ),
    plot = harmony_domains_by_sample,
    width = 16,
    height = 13,
    dpi = 300,
    bg = "white"
)

stopifnot(
    validObject(combined),
    ncol(Embeddings(combined, reduction = "pca")) == n_pcs_use,
    ncol(Embeddings(combined, reduction = "harmony")) == n_pcs_use,
    all(c("umap_pca", "umap_harmony") %in% Reductions(combined)),
    all(c("harmony_nn", "harmony_snn") %in% Graphs(combined)),
    "harmony_clusters_recomputed" %in% colnames(combined[[]]),
    is.factor(combined$harmony_clusters_recomputed),
    !anyNA(combined$harmony_clusters_recomputed),
    all(c("percent.mito", "percent.pltd") %in% colnames(combined[[]])),
    all(metadata_target_columns %in% colnames(combined[[]])),
    identical(as.character(combined$Barcode), colnames(combined)),
    is.factor(combined$sample_id),
    is.factor(combined$sample_domain),
    !anyNA(combined$sample_id),
    !anyNA(combined$sample_domain),
    identical(levels(combined$sample_id), sample_ids),
    identical(levels(combined$sample_domain), sample_domain_levels),
    all(vapply(
        factor_metadata_columns,
        function(x) is.factor(combined[[x, drop = TRUE]]),
        logical(1)
    ))
)

# Set the imported published Harmony cluster annotation as the active Seurat
# identity for exact reproduction of the downstream demonstration. The newly
# calculated clusters remain available as harmony_clusters_recomputed and should
# be used as the starting point when adapting the workflow to a new dataset.
Idents(combined) <- combined$harmony_clusters
stopifnot(identical(
    as.character(Idents(combined)),
    as.character(combined$harmony_clusters)
))

saveRDS(
    combined,
    file = file.path(
        "data",
        "processed",
        "maize_shoot_14samples_SCT_harmony_seurat_v5.rds"
    )
)
