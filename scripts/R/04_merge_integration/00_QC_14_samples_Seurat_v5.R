#!/usr/bin/env Rscript

# Quality-control summaries for the 14 maize shoot Visium datasets.
# Run from the maize_shoot_data_process_v2 repository root.
#
# Inputs:
#   data/processed/{sample_id}_seurat_v5.rds
#   data/metadata/metadata.csv
#   data/reference/maize_mitochondrial_genes.txt
#   data/reference/maize_plastid_genes.txt
#
# Outputs:
#   results/figures/04_QC/QC_spot_gene_histograms.png
#   results/figures/04_QC/QC_violin_by_sample.png
#   results/figures/04_QC/QC_features_vs_counts.png
#   results/figures/04_QC/QC_sample_domain_spot_counts.png
#   results/tables/QC_sample_summary.csv
#   results/tables/QC_gene_summary.csv
#   results/tables/QC_sample_domain_spot_counts.csv

suppressPackageStartupMessages({
    library(Seurat)
    library(Matrix)
    library(ggplot2)
    library(patchwork)
})

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

rds_files <- file.path(
    "data", "processed", paste0(sample_ids, "_seurat_v5.rds")
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

seurat_list <- lapply(seq_along(sample_ids), function(i) {
    object <- UpdateSeuratObject(readRDS(rds_files[i]))
    # SummarizedExperiment also exports Assays(), and it can mask the Seurat
    # method in an interactive session. Namespace this call so that it always
    # returns the character vector of Seurat assay names.
    assay_names <- SeuratObject::Assays(object)
    if (!"RNA" %in% assay_names && "Spatial" %in% assay_names) {
        object <- RenameAssays(object, Spatial = "RNA")
    }

    DefaultAssay(object) <- "RNA"
    object$sample_id <- sample_ids[i]
    object$orig.ident <- sample_ids[i]

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

combined <- merge(
    x = seurat_list[[1]],
    y = seurat_list[-1],
    add.cell.ids = sample_ids,
    project = "maize_shoot_14samples_QC",
    merge.data = FALSE,
    merge.dr = FALSE
)

# Rename merged cells to the format used by metadata.csv. The sample order in
# sample_ids defines the suffix: UL01 = _1_1 through DQ08 = _1_14. Preserve the
# original merged SAMPLE_BARCODE name in meta.data$old_colname.
old_cell_names <- colnames(combined)
cell_sample <- sub("_.*$", "", old_cell_names)
cell_barcode <- sub("^[^_]+_", "", old_cell_names)
cell_sample_number <- match(cell_sample, sample_ids)
stopifnot(!anyNA(cell_sample_number))

new_cell_names <- paste0(cell_barcode, "_1_", cell_sample_number)
if (anyDuplicated(new_cell_names)) {
    stop("The standardized cell names are not unique.")
}

combined$old_colname <- old_cell_names
combined <- RenameCells(combined, new.names = new_cell_names)

# Retain only non-overlapping shoot spots listed in metadata.csv, then align and
# add every metadata column by Barcode.
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
    "seurat_clusters", "harmony_clusters", "ms_ve", "domain_section","sample_domain",
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

rm(seurat_list)
gc()

figure_dir <- file.path("results", "figures", "04_QC")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path("results", "tables"), recursive = TRUE, showWarnings = FALSE)

spot_qc <- combined[[]]
spot_qc$spot_barcode <- rownames(spot_qc)

sample_summary <- do.call(rbind, lapply(sample_ids, function(sample_id) {
    x <- spot_qc[spot_qc$sample_id == sample_id, , drop = FALSE]
    data.frame(
        sample_id = sample_id,
        n_spots = nrow(x),
        median_detected_genes = median(x$nFeature_RNA),
        median_total_umi = median(x$nCount_RNA),
        median_percent_mito = median(x$percent.mito),
        median_percent_pltd = median(x$percent.pltd),
        max_percent_mito = max(x$percent.mito),
        max_percent_pltd = max(x$percent.pltd)
    )
}))
write.csv(
    sample_summary,
    file.path("results", "tables", "QC_sample_summary.csv"),
    row.names = FALSE
)

# Aggregate total reads and detected spots per gene without joining the 14
# sparse Seurat v5 count layers or creating another copy of the count matrix.
all_features <- rownames(combined[["RNA"]])
gene_total_umi <- setNames(numeric(length(all_features)), all_features)
gene_detected_spots <- setNames(numeric(length(all_features)), all_features)
count_layers <- Layers(combined[["RNA"]], search = "^counts")

for (layer_name in count_layers) {
    counts <- LayerData(combined, assay = "RNA", layer = layer_name)
    layer_features <- rownames(counts)
    gene_total_umi[layer_features] <- gene_total_umi[layer_features] +
        Matrix::rowSums(counts)
    gene_detected_spots[layer_features] <- gene_detected_spots[layer_features] +
        Matrix::rowSums(counts > 0)
}

gene_qc <- data.frame(
    gene = all_features,
    total_umi = unname(gene_total_umi[all_features]),
    detected_spots = unname(gene_detected_spots[all_features])
)
gene_qc$pass_minimum_100_reads <- gene_qc$total_umi >= 100
write.csv(
    gene_qc,
    file.path("results", "tables", "QC_gene_summary.csv"),
    row.names = FALSE
)

histogram_theme <- theme_gray(base_size = 11) +
    theme(plot.title = element_text(face = "bold"))
fill_colour <- "#F23846"

p_a <- ggplot(spot_qc, aes(x = nFeature_RNA)) +
    geom_histogram(bins = 45, fill = fill_colour, colour = NA) +
    labs(
        title = "Unique genes per tissue spot",
        x = "nFeature_RNA", y = "Count"
    ) +
    histogram_theme

p_b <- ggplot(spot_qc, aes(x = nCount_RNA)) +
    geom_histogram(bins = 45, fill = fill_colour, colour = NA) +
    labs(
        title = "Total mapped read counts per tissue spot",
        x = "nCount_RNA", y = "Count"
    ) +
    histogram_theme

p_c <- ggplot(gene_qc[gene_qc$total_umi > 0, ], aes(x = total_umi)) +
    geom_histogram(bins = 45, fill = fill_colour, colour = NA) +
    scale_x_log10() +
    labs(
        title = "Total mapped read counts per gene (log10 scale)",
        x = "Total UMI per gene", y = "Count"
    ) +
    histogram_theme

p_d <- ggplot(gene_qc, aes(x = detected_spots)) +
    geom_histogram(bins = 45, fill = fill_colour, colour = NA) +
    labs(
        title = "Total tissue spots per gene",
        x = "Detected tissue spots", y = "Count"
    ) +
    histogram_theme

qc_histograms <- (p_a | p_b) / (p_c | p_d) +
    plot_annotation(tag_levels = "A")
ggsave(
    file.path(figure_dir, "QC_spot_gene_histograms.png"),
    qc_histograms, width = 12, height = 10, dpi = 300, bg = "white"
)

qc_metrics <- c(
    nFeature_RNA = "Detected genes",
    nCount_RNA = "Total UMI counts",
    percent.mito = "Mitochondrial reads (%)",
    percent.pltd = "Plastid reads (%)"
)
qc_long <- do.call(rbind, lapply(names(qc_metrics), function(metric) {
    data.frame(
        sample_id = spot_qc$sample_id,
        metric = unname(qc_metrics[[metric]]),
        value = spot_qc[[metric]]
    )
}))
qc_long$metric <- factor(qc_long$metric, levels = unname(qc_metrics))

qc_violin <- ggplot(qc_long, aes(x = sample_id, y = value)) +
    geom_violin(fill = "grey92", colour = "grey45", scale = "width") +
    geom_boxplot(width = 0.12, outlier.shape = NA, linewidth = 0.25) +
    facet_wrap(~metric, scales = "free_y", ncol = 2) +
    labs(x = "Sample", y = NULL) +
    theme_classic(base_size = 10) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(
    file.path(figure_dir, "QC_violin_by_sample.png"),
    qc_violin, width = 12, height = 8, dpi = 300, bg = "white"
)

qc_scatter <- ggplot(
    spot_qc,
    aes(x = nCount_RNA, y = nFeature_RNA, colour = sample_id)
) +
    geom_point(size = 0.35, alpha = 0.35) +
    labs(
        title = "Detected genes versus total UMI counts",
        x = "Total UMI counts per tissue spot",
        y = "Detected genes per tissue spot",
        colour = "Sample"
    ) +
    theme_classic(base_size = 11)
ggsave(
    file.path(figure_dir, "QC_features_vs_counts.png"),
    qc_scatter, width = 8, height = 6, dpi = 300, bg = "white"
)

# Summarize the retained spots for every sample-domain combination. The
# complete factor levels retain zero-count combinations and keep the sample
# and anatomical-domain order stable across runs.
sample_domain_counts <- as.data.frame(
    table(
        sample_id = spot_qc$sample_id,
        domains = spot_qc$domains
    ),
    stringsAsFactors = FALSE,
    responseName = "n_spots"
)
sample_domain_counts$sample_id <- factor(
    sample_domain_counts$sample_id,
    levels = sample_ids
)
sample_domain_counts$domains <- factor(
    sample_domain_counts$domains,
    levels = allowed_domains
)
sample_domain_counts$sample_domain <- factor(
    paste(
        sample_domain_counts$sample_id,
        sample_domain_counts$domains,
        sep = "_"
    ),
    levels = sample_domain_levels
)
stopifnot(
    nrow(sample_domain_counts) == length(sample_domain_levels),
    !anyNA(sample_domain_counts$sample_domain),
    sum(sample_domain_counts$n_spots) == nrow(spot_qc)
)
write.csv(
    sample_domain_counts,
    file.path("results", "tables", "QC_sample_domain_spot_counts.csv"),
    row.names = FALSE
)

qc_sample_domain_heatmap <- ggplot(
    sample_domain_counts,
    aes(x = sample_id, y = domains, fill = n_spots)
) +
    geom_tile(colour = "white", linewidth = 0.4) +
    geom_text(
        aes(label = ifelse(n_spots > 0L, format(n_spots, big.mark = ","), "")),
        size = 3
    ) +
    scale_fill_gradient(
        low = "#F7FBFF",
        high = "#08519C",
        name = "Retained\nspots"
    ) +
    scale_y_discrete(limits = rev(allowed_domains)) +
    labs(
        title = "Retained tissue spots by sample and anatomical domain",
        x = "Sample",
        y = "Anatomical domain"
    ) +
    theme_classic(base_size = 11) +
    theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(face = "bold")
    )
ggsave(
    file.path(figure_dir, "QC_sample_domain_spot_counts.png"),
    qc_sample_domain_heatmap,
    width = 10,
    height = 5.5,
    dpi = 300,
    bg = "white"
)

message("QC complete: ", nrow(spot_qc), " spots and ", nrow(gene_qc), " genes.")
message("Genes with at least 100 total reads: ", sum(gene_qc$pass_minimum_100_reads))


session_dir <- file.path("results", "sessionInfo")
dir.create(session_dir, recursive = TRUE, showWarnings = FALSE)
writeLines(
    capture.output(sessionInfo()),
    file.path(session_dir, "04_merge_integration_sessionInfo.txt")
)
