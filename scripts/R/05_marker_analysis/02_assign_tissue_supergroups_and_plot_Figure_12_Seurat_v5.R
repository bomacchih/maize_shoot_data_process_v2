#!/usr/bin/env Rscript

# Assign the 33 Harmony clusters to anatomically defined tissue supergroups and
# reproduce Figure 12 (panels A-D) with Seurat v5.
#
# Run from the repository root:
#   source("scripts/R/05_marker_analysis/02_assign_tissue_supergroups_and_plot_Figure_12_Seurat_v5.R")
#
# Input:
#   data/processed/maize_shoot_14samples_SCT_harmony_seurat_v5.rds
#
# Optional upstream marker table:
#   results/tables/05_tissue_supergroups_Figure_12/markers_harmony_clusters_SCT_significant.csv
#   Generate it with 01_find_markers_harmony_clusters_SCT_Seurat_v5.R.
#
# Outputs:
#   results/tables/05_tissue_supergroups_Figure_12/cluster_to_tissue_supergroup.csv
#   results/tables/05_tissue_supergroups_Figure_12/cluster_domain_counts.csv
#   results/tables/05_tissue_supergroups_Figure_12/top10_markers_per_cluster_for_annotation.csv (if available)
#   results/figures/05_tissue_supergroups_Figure_12/Figure_12_A_clusters_UMAP.png
#   results/figures/05_tissue_supergroups_Figure_12/Figure_12_B_cluster_domain_heatmap.png
#   results/figures/05_tissue_supergroups_Figure_12/Figure_12_C_supergroups_UMAP.png
#   results/figures/05_tissue_supergroups_Figure_12/Figure_12_D_VR03_section2_spatial.png
#   results/figures/05_tissue_supergroups_Figure_12/Figure_12_composite.png
#   results/figures/05_tissue_supergroups_Figure_12/Figure_12_composite.pdf
#
# Assignment source:
# Supplementary Table 7-2, "The transfer of the unsupervised clusters to the
# supergroups of spatial information," in Supplementary_Tables_20251205.xlsx.
# The table defines 12 anatomical tissue supergroups plus sample_vari, which
# records sample-specific variation and is not treated as a tissue identity.
#
# DEMONSTRATION-ONLY WARNING:
# Harmony cluster numbers are analysis-specific labels, not transferable
# biological identities. The mapping below is valid only for the supplied maize
# shoot demonstration objects, for which the spot-level harmony_clusters values
# were verified to match exactly. For an independently processed dataset, users
# must identify their own clusters by examining marker genes, structural-domain
# distributions, spatial positions, and histology, and then replace the mapping
# below with a dataset-specific cluster-to-supergroup table.
# Step 04 also stores newly calculated clusters as harmony_clusters_recomputed.
# This Figure 12 reproduction intentionally does not use that field.

suppressPackageStartupMessages({
    library(Seurat)
    library(dplyr)
    library(ggplot2)
    library(patchwork)
})

input_file <- file.path(
    "data",
    "processed",
    "maize_shoot_14samples_SCT_harmony_seurat_v5.rds"
)
marker_file <- file.path(
    "results",
    "tables",
    "05_tissue_supergroups_Figure_12",
    "markers_harmony_clusters_SCT_significant.csv"
)
table_dir <- file.path(
    "results", "tables", "05_tissue_supergroups_Figure_12"
)
figure_dir <- file.path(
    "results", "figures", "05_tissue_supergroups_Figure_12"
)

dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

# UMAP display controls. The default values reproduce the orientation and
# coordinate range used in the manuscript figure. Change these values only if a
# newly calculated UMAP is mirrored or rotated relative to the reference figure.
umap_flip_x <- FALSE
umap_flip_y <- FALSE
umap_rotation_degrees <- 0
umap_x_limits <- c(-8.5, 11.0)
umap_y_limits <- c(-8.5, 8.5)

# Demonstration mapping for the supplied maize shoot dataset. Do not reuse these
# numeric cluster assignments for a newly integrated or reclustered dataset.
cluster_to_supergroup <- data.frame(
    cluster = as.character(0:32),
    spatial_pattern = c(
        "P5, most cell types except mature veins",
        "coleoptile, abaxial epidermis",
        "coleoptile mesophyll at the polar axis between vein and epidermis",
        "coleoptile cortex mesophyll at the equatorial axis",
        "coleoptile, abaxial epidermis",
        "SAM-P4, almost all cells",
        "coleoptile cortex mesophyll at the equatorial axis",
        "coleoptile cortex mesophyll at the equatorial axis",
        "SAM-P4, almost all cells except mature veins",
        "coleoptile cortex mesophyll at the equatorial axis, uneven",
        "P3-P5, mostly near mature veins",
        "coleoptile mesophyll next to the abaxial epidermis",
        "coleoptile cortex mesophyll at the equatorial axis, uneven",
        "coleoptile cortex mesophyll at the equatorial axis, uneven",
        "P1-P5, most cells, uneven",
        "coleoptile, adaxial epidermis",
        "coleoptile cortex mesophyll at the equatorial axis, uneven",
        "coleoptile cortex mesophyll surrounding veins",
        "P4-P5 mesophyll, uneven and sample-dependent",
        "coleoptile, abaxial epidermis",
        "abaxial side of the coleoptile vein (phloem)",
        "coleoptile mesophyll next to the abaxial epidermis",
        "P4-P5, most mesophyll, uneven",
        "P4, most mesophyll except veins",
        "mesophyll, sample-dependent",
        "adaxial side of the coleoptile vein (xylem)",
        "coleoptile cortex mesophyll at the equatorial axis, uneven",
        "coleoptile, abaxial epidermis",
        "not spatially specific",
        "not spatially specific",
        "not spatially specific",
        "not spatially specific",
        "not spatially specific"
    ),
    supergroup = c(
        "leaf_meso", "co_epi_ab", "co_meso_polar", "co_meso_core",
        "co_epi_ab", "SAM_P4", "co_meso_core", "co_meso_core",
        "leaf_meso", "co_meso_core", "leaf_vein", "co_meso_ab",
        "co_meso_core", "co_meso_core", "leaf_most", "co_epi_ad",
        "co_meso_core", "co_meso_nearvein", "sample_vari",
        "co_epi_ab", "co_vein_ab", "co_meso_ab", "leaf_most",
        "leaf_meso", "sample_vari", "co_vein_ad", "co_meso_core",
        "co_epi_ab", "sample_vari", "sample_vari", "sample_vari",
        "sample_vari", "sample_vari"
    ),
    stringsAsFactors = FALSE
)

supergroup_levels <- c(
    "SAM_P4",
    "leaf_meso",
    "leaf_vein",
    "co_meso_core",
    "co_meso_ab",
    "co_meso_polar",
    "co_meso_nearvein",
    "co_vein_ab",
    "co_vein_ad",
    "co_epi_ab",
    "co_epi_ad",
    "leaf_most",
    "sample_vari"
)

# Colors follow the visual identity of the manuscript figure. sample_vari is
# intentionally light gray because it is not an anatomical tissue supergroup.
supergroup_colors <- c(
    SAM_P4 = "#777777",
    leaf_meso = "#A8D98F",
    leaf_vein = "#E31A1C",
    co_meso_core = "#2C6DB2",
    co_meso_ab = "#79BFE2",
    co_meso_polar = "#25B7B2",
    co_meso_nearvein = "#2AA876",
    co_vein_ab = "#F29A8A",
    co_vein_ad = "#B78375",
    co_epi_ab = "#8E5AA9",
    co_epi_ad = "#D65AA5",
    leaf_most = "#C9E5B8",
    sample_vari = "#D0D0D0"
)

stopifnot(
    nrow(cluster_to_supergroup) == 33L,
    identical(cluster_to_supergroup$cluster, as.character(0:32)),
    setequal(unique(cluster_to_supergroup$supergroup), supergroup_levels),
    setequal(names(supergroup_colors), supergroup_levels)
)

combined <- readRDS(input_file)
required_metadata <- c(
    "sample_id", "section_id", "domains", "harmony_clusters"
)
stopifnot(
    inherits(combined, "Seurat"),
    "umap_harmony" %in% Reductions(combined),
    all(required_metadata %in% colnames(combined[[]])),
    "VR03" %in% Images(combined)
)

cluster_values <- trimws(as.character(combined$harmony_clusters))
unmapped_clusters <- setdiff(unique(cluster_values), cluster_to_supergroup$cluster)
if (length(unmapped_clusters) > 0L) {
    stop(
        "The following Harmony clusters are absent from Supplementary Table 7-2: ",
        paste(unmapped_clusters, collapse = ", ")
    )
}
missing_clusters <- setdiff(cluster_to_supergroup$cluster, unique(cluster_values))
if (length(missing_clusters) > 0L) {
    stop(
        "The following expected Harmony clusters are absent from the object: ",
        paste(missing_clusters, collapse = ", ")
    )
}

domain_levels <- c("SAM", "P1_P2", "P3", "P4", "P5", "coleoptile", "co_v")
domain_values <- trimws(as.character(combined$domains))
unexpected_domains <- setdiff(unique(domain_values), domain_levels)
if (length(unexpected_domains) > 0L || anyNA(domain_values)) {
    stop(
        "The domains metadata contains missing or unexpected values: ",
        paste(unexpected_domains, collapse = ", ")
    )
}

# Adding metadata does not change the active identity. Preserve and verify it.
active_ident_before <- Idents(combined)
supergroup_lookup <- setNames(
    cluster_to_supergroup$supergroup,
    cluster_to_supergroup$cluster
)
combined$supergroup <- factor(
    unname(supergroup_lookup[cluster_values]),
    levels = supergroup_levels
)
stopifnot(
    !anyNA(combined$supergroup),
    identical(
        unname(as.character(Idents(combined))),
        unname(as.character(active_ident_before))
    )
)

write.csv(
    cluster_to_supergroup,
    file.path(table_dir, "cluster_to_tissue_supergroup.csv"),
    row.names = FALSE
)

# Prepare manuscript-oriented UMAP coordinates. Both UMAP panels use the same
# transformed coordinates, ensuring that panels A and C are directly comparable.
umap_data <- as.data.frame(Embeddings(combined, reduction = "umap_harmony")[, 1:2])
colnames(umap_data) <- c("UMAP_1", "UMAP_2")
umap_data$spot <- rownames(umap_data)
umap_data$cluster <- factor(
    cluster_values[match(umap_data$spot, colnames(combined))],
    levels = as.character(0:32)
)
umap_data$supergroup <- factor(
    as.character(combined$supergroup)[match(umap_data$spot, colnames(combined))],
    levels = supergroup_levels
)

if (umap_flip_x) umap_data$UMAP_1 <- -umap_data$UMAP_1
if (umap_flip_y) umap_data$UMAP_2 <- -umap_data$UMAP_2

theta <- umap_rotation_degrees * pi / 180
umap_x_original <- umap_data$UMAP_1
umap_y_original <- umap_data$UMAP_2
umap_data$UMAP_1 <- (
    umap_x_original * cos(theta) - umap_y_original * sin(theta)
)
umap_data$UMAP_2 <- (
    umap_x_original * sin(theta) + umap_y_original * cos(theta)
)

cluster_colors <- setNames(
    grDevices::hcl.colors(33, palette = "Dark 3"),
    as.character(0:32)
)

base_umap_theme <- theme_classic(base_size = 10) +
    theme(
        axis.title = element_text(size = 11),
        axis.text = element_text(size = 8),
        legend.title = element_blank(),
        plot.margin = margin(5.5, 5.5, 5.5, 5.5)
    )

# Panel A: 33 unsupervised Harmony clusters.
pA <- ggplot(umap_data, aes(UMAP_1, UMAP_2, color = cluster)) +
    geom_point(size = 0.32, alpha = 0.95, stroke = 0) +
    scale_color_manual(values = cluster_colors, drop = FALSE) +
    coord_fixed(
        xlim = umap_x_limits,
        ylim = umap_y_limits,
        clip = "off"
    ) +
    labs(x = "UMAP-1", y = "UMAP-2") +
    guides(
        color = guide_legend(
            ncol = 2,
            override.aes = list(size = 2.2, alpha = 1)
        )
    ) +
    base_umap_theme +
    theme(legend.key.height = grid::unit(3.2, "mm"))

# Panel B: distribution of clusters across the seven structural domains.
cluster_domain_counts <- as.data.frame(
    table(
        cluster = factor(cluster_values, levels = as.character(0:32)),
        domain = factor(domain_values, levels = domain_levels)
    ),
    stringsAsFactors = FALSE
)
colnames(cluster_domain_counts)[3] <- "Count"
cluster_domain_counts$cluster <- factor(
    as.character(cluster_domain_counts$cluster),
    levels = as.character(0:32)
)
cluster_domain_counts$domain <- factor(
    as.character(cluster_domain_counts$domain),
    levels = domain_levels
)

write.csv(
    cluster_domain_counts,
    file.path(table_dir, "cluster_domain_counts.csv"),
    row.names = FALSE
)

pB <- ggplot(
    cluster_domain_counts,
    aes(x = domain, y = cluster, fill = Count)
) +
    geom_tile(color = "white", linewidth = 0.18) +
    scale_fill_gradient(low = "white", high = "#9E1B12") +
    scale_x_discrete(drop = FALSE) +
    scale_y_discrete(drop = FALSE) +
    labs(x = "Structural domain", y = "Harmony cluster", fill = "Count") +
    theme_classic(base_size = 9) +
    theme(
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        axis.ticks = element_blank(),
        legend.position = "right",
        plot.margin = margin(5.5, 5.5, 5.5, 5.5)
    )

# Panel C: the 12 anatomical tissue supergroups plus sample_vari.
label_positions <- umap_data %>%
    group_by(supergroup) %>%
    summarise(
        UMAP_1 = median(UMAP_1),
        UMAP_2 = median(UMAP_2),
        .groups = "drop"
    )

pC <- ggplot(umap_data, aes(UMAP_1, UMAP_2, color = supergroup)) +
    geom_point(size = 0.3, alpha = 0.9, stroke = 0) +
    scale_color_manual(values = supergroup_colors, drop = FALSE) +
    coord_fixed(
        xlim = umap_x_limits,
        ylim = umap_y_limits,
        clip = "off"
    ) +
    labs(x = "UMAP-1", y = "UMAP-2") +
    base_umap_theme +
    theme(legend.position = "none")

if (requireNamespace("ggrepel", quietly = TRUE)) {
    pC <- pC + ggrepel::geom_text_repel(
        data = label_positions,
        aes(label = supergroup),
        color = "black",
        size = 3.0,
        seed = 1234,
        box.padding = 0.35,
        point.padding = 0.15,
        min.segment.length = 0,
        max.overlaps = Inf,
        segment.color = "#555555",
        show.legend = FALSE
    )
} else {
    pC <- pC + geom_text(
        data = label_positions,
        aes(label = supergroup),
        color = "black",
        size = 2.8,
        check_overlap = TRUE,
        show.legend = FALSE
    )
}

# Panel D: spatial supergroup mapping for VR03 section 2. section_id is used
# instead of the generic Section2 label so only the intended biological sample
# is selected.
vr03_section2_cells <- rownames(combined[[]])[
    as.character(combined$section_id) == "VR03_S2"
]
if (length(vr03_section2_cells) == 0L) {
    stop("No spots were found for section_id == 'VR03_S2'.")
}

vr03_section2 <- subset(combined, cells = vr03_section2_cells)
stopifnot(
    all(as.character(vr03_section2$sample_id) == "VR03"),
    all(as.character(vr03_section2$section_id) == "VR03_S2")
)

pD <- SpatialDimPlot(
    object = vr03_section2,
    images = "VR03",
    group.by = "supergroup",
    cols = supergroup_colors,
    # The low-resolution image is the correctly scaled raster stored in the
    # combined Seurat object. Using "hires" misaligns the image and spots.
    image.scale = "lowres",
    image.alpha = 1,
    pt.size.factor = 3.5,
    alpha = c(1, 1),
    stroke = 0.1,
    crop = TRUE
) +
    theme(
        legend.position = "right",
        legend.title = element_text(size = 9),
        legend.text = element_text(size = 7),
        legend.key.height = grid::unit(3.5, "mm"),
        plot.margin = margin(5.5, 5.5, 5.5, 5.5)
    )

# If the marker-analysis script has already been run, make a compact table for
# manual comparison of marker genes, histology, and structural-domain labels.
if (file.exists(marker_file)) {
    markers <- read.csv(marker_file, check.names = FALSE)
    fold_change_column <- intersect(
        c("avg_log2FC", "avg_logFC"),
        colnames(markers)
    )
    fold_change_column <- fold_change_column[1]
    if (!is.na(fold_change_column) && "cluster" %in% colnames(markers)) {
        markers$cluster <- trimws(as.character(markers$cluster))
        markers$supergroup <- unname(supergroup_lookup[markers$cluster])
        top10_markers <- markers %>%
            arrange(
                cluster,
                desc(.data[[fold_change_column]]),
                p_val_adj
            ) %>%
            group_by(cluster) %>%
            slice_head(n = 10L) %>%
            ungroup()
        write.csv(
            top10_markers,
            file.path(
                table_dir,
                "top10_markers_per_cluster_for_annotation.csv"
            ),
            row.names = FALSE
        )
    } else {
        warning(
            "Marker table exists but lacks cluster or fold-change columns; ",
            "the top-marker annotation table was not generated."
        )
    }
} else {
    message(
        "Marker table not found. Run ",
        "01_find_markers_harmony_clusters_SCT_Seurat_v5.R before manual ",
        "marker review."
    )
}

# Save individual panels for flexible manuscript layout.
ggsave(
    file.path(figure_dir, "Figure_12_A_clusters_UMAP.png"),
    pA,
    width = 7.2,
    height = 6.3,
    dpi = 600,
    bg = "white"
)
ggsave(
    file.path(figure_dir, "Figure_12_B_cluster_domain_heatmap.png"),
    pB,
    width = 4.2,
    height = 6.3,
    dpi = 600,
    bg = "white"
)
ggsave(
    file.path(figure_dir, "Figure_12_C_supergroups_UMAP.png"),
    pC,
    width = 7.2,
    height = 6.3,
    dpi = 600,
    bg = "white"
)
ggsave(
    file.path(figure_dir, "Figure_12_D_VR03_section2_spatial.png"),
    pD,
    width = 7.0,
    height = 6.3,
    dpi = 600,
    bg = "white"
)

figure_12 <- (
    (pA + pB + plot_layout(widths = c(1.7, 1.0))) /
        (pC + pD + plot_layout(widths = c(1.15, 1.0)))
) +
    plot_annotation(tag_levels = "A") &
    theme(plot.tag = element_text(face = "bold", size = 16))

ggsave(
    file.path(figure_dir, "Figure_12_composite.png"),
    figure_12,
    width = 14.5,
    height = 12.0,
    dpi = 600,
    bg = "white"
)
ggsave(
    file.path(figure_dir, "Figure_12_composite.pdf"),
    figure_12,
    width = 14.5,
    height = 12.0,
    device = grDevices::pdf,
    bg = "white"
)

message(
    "Figure 12 workflow complete. Assigned ",
    ncol(combined),
    " spots from 33 clusters to 12 anatomical tissue supergroups plus ",
    "sample_vari. Panel D contains ",
    length(vr03_section2_cells),
    " spots from VR03 section 2."
)


session_dir <- file.path("results", "sessionInfo")
dir.create(session_dir, recursive = TRUE, showWarnings = FALSE)
writeLines(
    capture.output(sessionInfo()),
    file.path(session_dir, "05_marker_analysis_sessionInfo.txt")
)
