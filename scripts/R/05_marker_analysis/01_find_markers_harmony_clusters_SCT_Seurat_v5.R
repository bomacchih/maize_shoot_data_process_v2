#!/usr/bin/env Rscript

# Identify positive marker genes for Harmony clusters in the final maize shoot
# Seurat v5 object.
#
# Run this script from the repository root:
#   source("scripts/R/05_marker_analysis/01_find_markers_harmony_clusters_SCT_Seurat_v5.R")
#
# Input:
#   data/processed/maize_shoot_14samples_SCT_harmony_seurat_v5.rds
#
# Outputs:
#   results/tables/05_tissue_supergroups_Figure_12/markers_harmony_clusters_SCT_all.csv
#   results/tables/05_tissue_supergroups_Figure_12/markers_harmony_clusters_SCT_significant.csv
#   results/tables/05_tissue_supergroups_Figure_12/marker_genes_harmony_clusters_SCT_significant.txt
#   results/logs/markers_harmony_clusters_SCT_summary.txt
#
# This is an exploratory spot-level marker screen for cluster characterization,
# matching the original analysis settings. It is not the replicate-level
# differential-expression analysis described in the revised manuscript. Spots
# from the same section or capture area are not independent biological
# replicates; use raw-count pseudobulk methods for inferential domain comparisons.
#
# DEMONSTRATION-ONLY IDENTITY:
# This script intentionally uses harmony_clusters, the published cluster labels
# imported from metadata.csv, to reproduce the manuscript marker analysis and
# Figure 12 mapping. Step 04 also calculates harmony_clusters_recomputed from the
# current Harmony graph. Users analyzing their own data should characterize the
# recomputed clusters and must not transfer the published numeric labels.

suppressPackageStartupMessages({
    library(Seurat)
    library(dplyr)
})

input_file <- file.path(
    "data",
    "processed",
    "maize_shoot_14samples_SCT_harmony_seurat_v5.rds"
)
table_dir <- file.path(
    "results", "tables", "05_tissue_supergroups_Figure_12"
)
log_dir <- file.path("results", "logs")

# Marker-screen parameters retained from the original analysis script:
# - min.pct = 0.25 requires detection in at least 25% of spots in either group.
# - logfc.threshold = 0.25 removes very small positive effects before testing.
# - adjusted p-value < 0.05 defines the retained exploratory marker table.
# These values must not be confused with the manuscript's pseudobulk DEG
# criteria (absolute log2 fold change > 1 and FDR < 0.05).
marker_min_pct <- 0.25
marker_logfc_threshold <- 0.25
marker_adjusted_p_cutoff <- 0.05

dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

combined <- readRDS(input_file)
# SummarizedExperiment also exports Assays() and may mask Seurat's accessor in
# an interactive session, so explicitly use the SeuratObject implementation.
assay_names <- SeuratObject::Assays(combined)
stopifnot(
    inherits(combined, "Seurat"),
    "SCT" %in% assay_names,
    "harmony_clusters" %in% colnames(combined[[]]),
    is.factor(combined$harmony_clusters)
)

# Use the imported Harmony cluster annotation as the active identity. The
# general object name keeps this script reusable and avoids hard-coded object
# names from the original analysis workspace.
Idents(combined) <- combined$harmony_clusters
stopifnot(identical(
    unname(as.character(Idents(combined))),
    unname(as.character(combined$harmony_clusters))
))

DefaultAssay(combined) <- "SCT"

# The SCT assay contains one model for each sample. Recalculate comparable SCT
# counts using the minimum median UMI before exploratory marker testing.
combined <- PrepSCTFindMarkers(
    object = combined,
    assay = "SCT",
    verbose = TRUE
)

set.seed(1234)
markers_SCT_all <- FindAllMarkers(
    object = combined,
    assay = "SCT",
    only.pos = TRUE,
    min.pct = marker_min_pct,
    logfc.threshold = marker_logfc_threshold,
    test.use = "wilcox",
    return.thresh = 1,
    verbose = TRUE
)

if (nrow(markers_SCT_all) == 0L) {
    stop("FindAllMarkers() returned no marker genes.")
}

markers_SCT_all <- markers_SCT_all %>%
    arrange(cluster, p_val_adj, p_val)

significant_markers_SCT_all <- markers_SCT_all %>%
    filter(
        !is.na(p_val_adj),
        p_val_adj < marker_adjusted_p_cutoff
    )

significant_marker_genes_SCT_all <- significant_markers_SCT_all %>%
    pull(gene) %>%
    unique()

write.csv(
    markers_SCT_all,
    file = file.path(table_dir, "markers_harmony_clusters_SCT_all.csv"),
    row.names = FALSE
)
write.csv(
    significant_markers_SCT_all,
    file = file.path(
        table_dir,
        "markers_harmony_clusters_SCT_significant.csv"
    ),
    row.names = FALSE
)
writeLines(
    significant_marker_genes_SCT_all,
    con = file.path(
        table_dir,
        "marker_genes_harmony_clusters_SCT_significant.txt"
    )
)

significant_by_cluster <- significant_markers_SCT_all %>%
    count(cluster, name = "significant_markers")

summary_lines <- c(
    paste0("Input object: ", input_file),
    paste0("Spots: ", ncol(combined)),
    paste0("Features: ", nrow(combined)),
    paste0("Identity field: harmony_clusters"),
    paste0("Identity source: published labels imported from metadata.csv"),
    paste0("Identity levels: ", paste(levels(Idents(combined)), collapse = ", ")),
    "Assay: SCT",
    "Test: Wilcoxon rank-sum",
    "only.pos: TRUE",
    paste0("min.pct: ", marker_min_pct),
    paste0("logfc.threshold: ", marker_logfc_threshold),
    paste0(
        "Exploratory marker criterion: adjusted p-value < ",
        marker_adjusted_p_cutoff
    ),
    paste0(
        "Separate pseudobulk DEG criterion: absolute log2 fold change > 1 ",
        "and FDR < 0.05"
    ),
    paste0("Candidate marker rows: ", nrow(markers_SCT_all)),
    paste0("Significant marker rows: ", nrow(significant_markers_SCT_all)),
    paste0("Unique significant genes: ", length(significant_marker_genes_SCT_all)),
    "",
    "Significant marker rows by cluster:",
    paste0(
        significant_by_cluster$cluster,
        ": ",
        significant_by_cluster$significant_markers
    )
)
writeLines(
    summary_lines,
    con = file.path(log_dir, "markers_harmony_clusters_SCT_summary.txt")
)

message(
    "Marker analysis complete: ", nrow(significant_markers_SCT_all),
    " significant marker rows representing ",
    length(significant_marker_genes_SCT_all), " unique genes."
)


session_dir <- file.path("results", "sessionInfo")
dir.create(session_dir, recursive = TRUE, showWarnings = FALSE)
writeLines(
    capture.output(sessionInfo()),
    file.path(session_dir, "05_marker_analysis_sessionInfo.txt")
)
