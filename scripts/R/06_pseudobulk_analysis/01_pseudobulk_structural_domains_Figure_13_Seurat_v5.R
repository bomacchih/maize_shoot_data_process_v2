#!/usr/bin/env Rscript

# Generate replicate-level pseudobulk profiles for the seven maize shoot
# structural domains and reproduce Figure 13 (domain-mean PCA and hierarchical
# clustering) with Seurat v5 and edgeR.
#
# Run from the repository root:
#   source("scripts/R/06_pseudobulk_analysis/01_pseudobulk_structural_domains_Figure_13_Seurat_v5.R")
#
# Experimental unit:
#   sample_id is the biological-replicate identifier. Raw UMI counts from all
#   spots and serial sections belonging to the same sample_id and domain are
#   summed into one pseudobulk library. Individual spots and sections are not
#   treated as independent biological replicates.
#
# Input:
#   data/processed/maize_shoot_14samples_SCT_harmony_seurat_v5.rds
#
# Main outputs:
#   results/tables/Figure_13/pseudobulk_library_metadata.csv
#   results/tables/Figure_13/pseudobulk_raw_counts.csv.gz
#   results/tables/Figure_13/pseudobulk_TMM_log2CPM.csv.gz
#   results/tables/Figure_13/domain_mean_TMM_log2CPM.csv.gz
#   results/tables/Figure_13/replicate_PCA_coordinates.csv
#   results/tables/Figure_13/domain_mean_PCA_coordinates.csv
#   results/figures/Figure_13/Figure_13_A_domain_mean_PCA.png
#   results/figures/Figure_13/Figure_13_B_domain_mean_hierarchical_clustering.png
#   results/figures/Figure_13/Figure_13_composite.png
#   results/figures/Figure_13/Figure_13_composite.pdf
#   results/figures/Figure_13/pseudobulk_replicate_PCA_diagnostic.png
#   results/sessionInfo/06_pseudobulk_analysis_sessionInfo.txt

suppressPackageStartupMessages({
    library(Seurat)
    library(SeuratObject)
    library(Matrix)
    library(edgeR)
    library(dplyr)
    library(ggplot2)
    library(patchwork)
})

set.seed(1234)

input_file <- file.path(
    "data",
    "processed",
    "maize_shoot_14samples_SCT_harmony_seurat_v5.rds"
)
table_dir <- file.path("results", "tables", "Figure_13")
figure_dir <- file.path("results", "figures", "Figure_13")
log_dir <- file.path("results", "logs")
session_dir <- file.path("results", "sessionInfo")

dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(session_dir, recursive = TRUE, showWarnings = FALSE)

# These metadata columns define the experimental unit and anatomical group.
replicate_column <- "sample_id"
domain_column <- "domains"
domain_levels <- c(
    "SAM", "P1_P2", "P3", "P4", "P5", "coleoptile", "co_v"
)

domain_labels <- c(
    SAM = "SAM",
    P1_P2 = "P1-P2",
    P3 = "P3",
    P4 = "P4",
    P5 = "P5",
    coleoptile = "coleoptile",
    co_v = "co-v"
)

domain_colors <- c(
    SAM = "#F15A4F",
    P1_P2 = "#C51B3A",
    P3 = "#14863D",
    P4 = "#18B6B2",
    P5 = "#2F68B2",
    coleoptile = "#65539B",
    co_v = "#5C2B6E"
)

# Optional scRNA-seq comparison. Set this to a CSV path to enable the
# correlation step. The first column must contain gene identifiers, and the
# remaining columns must contain normalized reference profiles.
scrna_reference_file <- NULL

write_matrix_csv_gz <- function(matrix_to_write, output_file) {
    output_connection <- gzfile(output_file, open = "wt")
    on.exit(close(output_connection), add = TRUE)

    output_table <- data.frame(
        gene = rownames(matrix_to_write),
        as.data.frame(as.matrix(matrix_to_write), check.names = FALSE),
        check.names = FALSE
    )
    utils::write.csv(
        output_table,
        file = output_connection,
        row.names = FALSE,
        quote = FALSE
    )
}

calculate_spearman_correlations <- function(
    spatial_profiles,
    scrna_profiles,
    minimum_common_genes = 100L
) {
    common_genes <- intersect(
        rownames(spatial_profiles),
        rownames(scrna_profiles)
    )

    if (length(common_genes) < minimum_common_genes) {
        stop(
            "Only ",
            length(common_genes),
            " genes overlap between the spatial and scRNA-seq profiles; ",
            "at least ",
            minimum_common_genes,
            " are required."
        )
    }

    stats::cor(
        x = spatial_profiles[common_genes, , drop = FALSE],
        y = scrna_profiles[common_genes, , drop = FALSE],
        method = "spearman",
        use = "pairwise.complete.obs"
    )
}

# Convert an hclust object to line segments so the dendrogram can be drawn
# using ggplot2 without an additional plotting package.
hclust_to_ggplot_data <- function(hclust_object) {
    dendrogram_object <- as.dendrogram(hclust_object)
    leaf_labels <- labels(dendrogram_object)
    leaf_x <- stats::setNames(seq_along(leaf_labels), leaf_labels)
    segment_table <- data.frame(
        x = numeric(),
        y = numeric(),
        xend = numeric(),
        yend = numeric()
    )

    walk_dendrogram <- function(node) {
        node_height <- attr(node, "height")

        if (is.leaf(node)) {
            node_label <- attr(node, "label")
            return(list(x = unname(leaf_x[[node_label]]), height = 0))
        }

        child_positions <- lapply(node, walk_dendrogram)
        child_x <- vapply(child_positions, `[[`, numeric(1), "x")
        child_height <- vapply(
            child_positions,
            `[[`,
            numeric(1),
            "height"
        )
        node_x <- mean(child_x)

        for (child_index in seq_along(child_x)) {
            segment_table <<- rbind(
                segment_table,
                data.frame(
                    x = child_x[child_index],
                    y = child_height[child_index],
                    xend = child_x[child_index],
                    yend = node_height
                )
            )
        }

        segment_table <<- rbind(
            segment_table,
            data.frame(
                x = min(child_x),
                y = node_height,
                xend = max(child_x),
                yend = node_height
            )
        )

        list(x = node_x, height = node_height)
    }

    walk_dendrogram(dendrogram_object)
    list(segments = segment_table, labels = leaf_labels)
}

make_leaf_ellipse <- function(pca_coordinates, leaf_domains, scale = 1.7) {
    leaf_points <- pca_coordinates[
        pca_coordinates$domain %in% leaf_domains,
        c("PC1", "PC2"),
        drop = FALSE
    ]

    if (nrow(leaf_points) < 3L) {
        return(NULL)
    }

    covariance_matrix <- stats::cov(leaf_points)
    eigen_result <- eigen(covariance_matrix, symmetric = TRUE)
    eigen_values <- pmax(eigen_result$values, 0)
    angles <- seq(0, 2 * pi, length.out = 201L)
    unit_circle <- cbind(cos(angles), sin(angles))
    ellipse_points <- unit_circle %*%
        diag(sqrt(eigen_values), nrow = 2L) %*%
        t(eigen_result$vectors)
    ellipse_points <- ellipse_points * scale
    ellipse_points <- sweep(
        ellipse_points,
        MARGIN = 2L,
        STATS = colMeans(leaf_points),
        FUN = "+"
    )

    data.frame(PC1 = ellipse_points[, 1], PC2 = ellipse_points[, 2])
}

if (!file.exists(input_file)) {
    stop("Input Seurat object not found: ", input_file)
}

combined <- readRDS(input_file)
active_ident_before <- Idents(combined)

required_metadata <- c(replicate_column, domain_column)
missing_metadata <- setdiff(required_metadata, colnames(combined[[]]))
if (length(missing_metadata) > 0L) {
    stop(
        "Required metadata columns are missing: ",
        paste(missing_metadata, collapse = ", ")
    )
}
if (!"RNA" %in% Assays(combined)) {
    stop("The input object does not contain an RNA assay.")
}

# Join only the raw count layers. This creates one counts layer while
# preserving the original raw UMI values.
count_layers <- Layers(combined[["RNA"]], search = "^counts")
if (length(count_layers) == 0L) {
    stop("No raw RNA count layers were found.")
}
if (length(count_layers) > 1L) {
    combined <- JoinLayers(
        object = combined,
        assay = "RNA"
    )
}

if (!"counts" %in% Layers(combined[["RNA"]])) {
    stop(
        "Joining the RNA assay did not produce a single raw 'counts' layer."
    )
}

raw_counts <- LayerData(
    object = combined,
    assay = "RNA",
    layer = "counts"
)
if (!inherits(raw_counts, "sparseMatrix")) {
    raw_counts <- as(raw_counts, "dgCMatrix")
}

# Align metadata to the raw-count columns and retain only the seven domains.
spot_metadata <- combined[[]][colnames(raw_counts), , drop = FALSE]
spot_metadata$biological_replicate <- trimws(
    as.character(spot_metadata[[replicate_column]])
)
spot_metadata$domain <- trimws(
    as.character(spot_metadata[[domain_column]])
)

valid_spots <- !is.na(spot_metadata$biological_replicate) &
    nzchar(spot_metadata$biological_replicate) &
    spot_metadata$domain %in% domain_levels

raw_counts <- raw_counts[, valid_spots, drop = FALSE]
spot_metadata <- spot_metadata[valid_spots, , drop = FALSE]

if (ncol(raw_counts) == 0L) {
    stop("No spots remained after validating replicate and domain metadata.")
}
if (!identical(colnames(raw_counts), rownames(spot_metadata))) {
    stop("Raw count columns and spot metadata rows are not aligned.")
}

# Define one pseudobulk library per biological replicate and structural domain.
spot_metadata$group_id <- paste(
    spot_metadata$biological_replicate,
    spot_metadata$domain,
    sep = "__"
)

group_metadata <- unique(spot_metadata[, c(
    "group_id",
    "biological_replicate",
    "domain"
)])
replicate_levels <- unique(spot_metadata$biological_replicate)
group_metadata <- group_metadata[order(
    match(group_metadata$biological_replicate, replicate_levels),
    match(group_metadata$domain, domain_levels)
), , drop = FALSE]
rownames(group_metadata) <- group_metadata$group_id

group_factor <- factor(
    spot_metadata$group_id,
    levels = group_metadata$group_id
)

# Sparse spot-by-library membership matrix; multiplying the gene-by-spot raw
# count matrix by this matrix sums raw UMIs without creating a dense spot matrix.
membership_matrix <- sparseMatrix(
    i = seq_along(group_factor),
    j = as.integer(group_factor),
    x = 1,
    dims = c(length(group_factor), nlevels(group_factor)),
    dimnames = list(colnames(raw_counts), levels(group_factor))
)

pseudobulk_counts <- as.matrix(raw_counts %*% membership_matrix)
storage.mode(pseudobulk_counts) <- "numeric"
group_metadata$spot_count <- as.integer(table(group_factor)[
    group_metadata$group_id
])

if (any(colSums(pseudobulk_counts) == 0)) {
    stop("At least one pseudobulk library has a total raw count of zero.")
}

# TMM normalization is calculated across replicate-domain pseudobulk libraries.
# edgeR::cpm(..., log = TRUE) returns log2 counts per million.
dge <- DGEList(
    counts = pseudobulk_counts,
    samples = group_metadata
)
if ("normLibSizes" %in% getNamespaceExports("edgeR")) {
    dge <- normLibSizes(dge, method = "TMM")
} else {
    # Compatibility with edgeR releases before normLibSizes() was introduced.
    dge <- calcNormFactors(dge, method = "TMM")
}
pseudobulk_log2cpm <- cpm(dge, log = TRUE, prior.count = 1)

group_metadata$library_size <- dge$samples$lib.size
group_metadata$TMM_normalization_factor <- dge$samples$norm.factors
group_metadata$effective_library_size <- with(
    group_metadata,
    library_size * TMM_normalization_factor
)

# The mean gives equal weight to each biological replicate represented in a
# domain. These domain means are used only for descriptive PCA and clustering.
missing_domains <- setdiff(domain_levels, unique(group_metadata$domain))
if (length(missing_domains) > 0L) {
    stop(
        "No pseudobulk profiles were available for: ",
        paste(missing_domains, collapse = ", ")
    )
}

domain_mean_log2cpm <- vapply(
    domain_levels,
    function(current_domain) {
        domain_columns <- group_metadata$group_id[
            group_metadata$domain == current_domain
        ]
        rowMeans(
            pseudobulk_log2cpm[, domain_columns, drop = FALSE]
        )
    },
    numeric(nrow(pseudobulk_log2cpm))
)
rownames(domain_mean_log2cpm) <- rownames(pseudobulk_log2cpm)
colnames(domain_mean_log2cpm) <- domain_levels

# Export the profiles and library metadata before dimensional reduction.
write.csv(
    group_metadata,
    file.path(table_dir, "pseudobulk_library_metadata.csv"),
    row.names = FALSE
)
write_matrix_csv_gz(
    pseudobulk_counts,
    file.path(table_dir, "pseudobulk_raw_counts.csv.gz")
)
write_matrix_csv_gz(
    pseudobulk_log2cpm,
    file.path(table_dir, "pseudobulk_TMM_log2CPM.csv.gz")
)
write_matrix_csv_gz(
    domain_mean_log2cpm,
    file.path(table_dir, "domain_mean_TMM_log2CPM.csv.gz")
)

# Replicate-level PCA is a diagnostic for within-domain replicate variation.
replicate_pca <- prcomp(
    t(pseudobulk_log2cpm),
    center = TRUE,
    scale. = FALSE
)
replicate_variance <- 100 * replicate_pca$sdev^2 /
    sum(replicate_pca$sdev^2)
replicate_pca_coordinates <- data.frame(
    group_metadata,
    PC1 = replicate_pca$x[, 1],
    PC2 = replicate_pca$x[, 2],
    row.names = NULL,
    check.names = FALSE
)
replicate_pca_coordinates$domain <- factor(
    replicate_pca_coordinates$domain,
    levels = domain_levels
)

write.csv(
    replicate_pca_coordinates,
    file.path(table_dir, "replicate_PCA_coordinates.csv"),
    row.names = FALSE
)

p_replicates <- ggplot(
    replicate_pca_coordinates,
    aes(PC1, PC2, color = domain, label = biological_replicate)
) +
    geom_point(size = 2.2, alpha = 0.85) +
    scale_color_manual(
        values = domain_colors,
        breaks = domain_levels,
        labels = domain_labels,
        drop = FALSE
    ) +
    labs(
        x = sprintf("PC1 (%.1f%% variance explained)", replicate_variance[1]),
        y = sprintf("PC2 (%.1f%% variance explained)", replicate_variance[2]),
        color = "Domain",
        title = "Replicate-level pseudobulk profiles"
    ) +
    theme_classic(base_size = 11) +
    theme(legend.position = "right")

if (requireNamespace("ggrepel", quietly = TRUE)) {
    p_replicates <- p_replicates + ggrepel::geom_text_repel(
        size = 2.5,
        seed = 1234,
        max.overlaps = Inf,
        show.legend = FALSE
    )
}

ggsave(
    file.path(figure_dir, "pseudobulk_replicate_PCA_diagnostic.png"),
    p_replicates,
    width = 8.2,
    height = 6.4,
    dpi = 600,
    bg = "white"
)

# Figure 13A: descriptive PCA of the seven equal-weight domain mean profiles.
domain_pca <- prcomp(
    t(domain_mean_log2cpm),
    center = TRUE,
    scale. = FALSE
)
domain_variance <- 100 * domain_pca$sdev^2 / sum(domain_pca$sdev^2)
domain_pca_coordinates <- data.frame(
    domain = factor(domain_levels, levels = domain_levels),
    PC1 = domain_pca$x[domain_levels, 1],
    PC2 = domain_pca$x[domain_levels, 2],
    label = unname(domain_labels[domain_levels]),
    row.names = NULL
)

write.csv(
    domain_pca_coordinates,
    file.path(table_dir, "domain_mean_PCA_coordinates.csv"),
    row.names = FALSE
)

leaf_domains <- c("P1_P2", "P3", "P4", "P5")
leaf_ellipse <- make_leaf_ellipse(
    domain_pca_coordinates,
    leaf_domains = leaf_domains
)
sam_coordinates <- domain_pca_coordinates[
    domain_pca_coordinates$domain == "SAM",
    c("PC1", "PC2")
]
leaf_centroid <- colMeans(domain_pca_coordinates[
    domain_pca_coordinates$domain %in% leaf_domains,
    c("PC1", "PC2")
])

pA <- ggplot(
    domain_pca_coordinates,
    aes(PC1, PC2, color = domain)
) +
    geom_point(size = 4.0) +
    scale_color_manual(
        values = domain_colors,
        breaks = domain_levels,
        labels = domain_labels,
        drop = FALSE
    ) +
    labs(
        x = sprintf("PC1 (%.1f%% variance explained)", domain_variance[1]),
        y = sprintf("PC2 (%.1f%% variance explained)", domain_variance[2]),
        color = NULL
    ) +
    theme_classic(base_size = 11) +
    theme(
        legend.position = "right",
        legend.key.height = grid::unit(4, "mm")
    )

if (!is.null(leaf_ellipse)) {
    pA <- pA + geom_path(
        data = leaf_ellipse,
        aes(PC1, PC2),
        inherit.aes = FALSE,
        color = "#2C75B8",
        linewidth = 0.6,
        linetype = "dashed"
    )
}

pA <- pA +
    geom_curve(
        data = data.frame(
            x = sam_coordinates$PC1,
            y = sam_coordinates$PC2,
            xend = leaf_centroid[["PC1"]],
            yend = leaf_centroid[["PC2"]]
        ),
        aes(x = x, y = y, xend = xend, yend = yend),
        inherit.aes = FALSE,
        color = "grey70",
        linewidth = 2.4,
        curvature = -0.18,
        arrow = grid::arrow(
            length = grid::unit(4, "mm"),
            type = "closed"
        )
    )

if (requireNamespace("ggrepel", quietly = TRUE)) {
    pA <- pA + ggrepel::geom_text_repel(
        aes(label = label),
        color = "black",
        size = 3.2,
        seed = 1234,
        max.overlaps = Inf,
        show.legend = FALSE
    )
} else {
    pA <- pA + geom_text(
        aes(label = label),
        color = "black",
        size = 3.0,
        nudge_y = 0.08 * diff(range(domain_pca_coordinates$PC2)),
        check_overlap = TRUE,
        show.legend = FALSE
    )
}

# Figure 13B: Euclidean distance and complete-linkage clustering of the same
# seven domain mean profiles used in panel A.
domain_distance <- dist(t(domain_mean_log2cpm), method = "euclidean")
domain_hclust <- hclust(domain_distance, method = "complete")
dendrogram_data <- hclust_to_ggplot_data(domain_hclust)
dendrogram_label_table <- data.frame(
    x = seq_along(dendrogram_data$labels),
    label = unname(domain_labels[dendrogram_data$labels])
)

pB <- ggplot(dendrogram_data$segments) +
    geom_segment(
        aes(x = x, y = y, xend = xend, yend = yend),
        linewidth = 0.8,
        lineend = "square"
    ) +
    scale_x_continuous(
        breaks = dendrogram_label_table$x,
        labels = dendrogram_label_table$label,
        expand = expansion(mult = c(0.08, 0.08))
    ) +
    labs(x = NULL, y = "Height") +
    theme_classic(base_size = 11) +
    theme(
        axis.text.x = element_text(face = "bold", angle = 0, vjust = 1),
        axis.ticks.x = element_blank(),
        axis.line.x = element_blank()
    )

ggsave(
    file.path(figure_dir, "Figure_13_A_domain_mean_PCA.png"),
    pA,
    width = 7.2,
    height = 6.2,
    dpi = 600,
    bg = "white"
)
ggsave(
    file.path(
        figure_dir,
        "Figure_13_B_domain_mean_hierarchical_clustering.png"
    ),
    pB,
    width = 6.0,
    height = 6.2,
    dpi = 600,
    bg = "white"
)

figure_13 <- (pA + pB) +
    plot_layout(widths = c(1.15, 1.0)) +
    plot_annotation(tag_levels = "A") &
    theme(plot.tag = element_text(face = "bold", size = 16))

ggsave(
    file.path(figure_dir, "Figure_13_composite.png"),
    figure_13,
    width = 13.2,
    height = 6.2,
    dpi = 600,
    bg = "white"
)
ggsave(
    file.path(figure_dir, "Figure_13_composite.pdf"),
    figure_13,
    width = 13.2,
    height = 6.2,
    device = grDevices::pdf,
    bg = "white"
)

# Optional comparison with normalized scRNA-seq reference profiles.
if (!is.null(scrna_reference_file)) {
    if (!file.exists(scrna_reference_file)) {
        stop("scRNA-seq reference file not found: ", scrna_reference_file)
    }

    scrna_table <- read.csv(
        scrna_reference_file,
        check.names = FALSE,
        stringsAsFactors = FALSE
    )
    if (ncol(scrna_table) < 2L) {
        stop(
            "The scRNA-seq reference CSV must contain a gene column and at ",
            "least one normalized profile column."
        )
    }

    scrna_profiles <- as.matrix(scrna_table[, -1, drop = FALSE])
    rownames(scrna_profiles) <- as.character(scrna_table[[1]])
    storage.mode(scrna_profiles) <- "numeric"

    spearman_correlations <- calculate_spearman_correlations(
        spatial_profiles = domain_mean_log2cpm,
        scrna_profiles = scrna_profiles
    )
    write.csv(
        spearman_correlations,
        file.path(
            table_dir,
            "spatial_domain_vs_scRNAseq_Spearman_correlations.csv"
        ),
        row.names = TRUE
    )
}

stopifnot(
    identical(
        unname(as.character(Idents(combined))),
        unname(as.character(active_ident_before))
    )
)

session_information <- capture.output(sessionInfo())
writeLines(
    session_information,
    file.path(session_dir, "06_pseudobulk_analysis_sessionInfo.txt")
)

message(
    "Figure 13 pseudobulk workflow complete: ",
    ncol(pseudobulk_counts),
    " replicate-domain libraries from ",
    length(unique(group_metadata$biological_replicate)),
    " biological replicates and seven structural domains."
)
