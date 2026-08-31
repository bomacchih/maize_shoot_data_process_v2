#!/usr/bin/env Rscript

# Postprocess the original single-time-series TO-GCN output, generate the
# TF-level heatmap and network/profile panels, and export level-specific gene
# sets for AgriGO or optional local Fisher-exact enrichment.
#
# Run after the R preprocessing and shell TO-GCN steps:
#   source("scripts/R/11_TO_GCN/02_postprocess_TO_GCN_plot_and_GO.R")

suppressPackageStartupMessages({
    library(dplyr)
    library(tidyr)
    library(ggplot2)
    library(patchwork)
})

set.seed(1)

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

read_expression_matrix <- function(path) {
    x <- read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
    gene_column <- detect_column(names(x), c("gene_id", "gene"))
    if (is.na(gene_column)) gene_column <- names(x)[1L]
    gene_ids <- trimws(as.character(x[[gene_column]]))
    x[[gene_column]] <- NULL
    matrix_out <- as.matrix(x)
    storage.mode(matrix_out) <- "numeric"
    rownames(matrix_out) <- gene_ids
    matrix_out
}

clamp_values <- function(x, range) {
    pmax(range[1L], pmin(range[2L], x))
}

project_root <- find_project_root()
input_dir <- file.path(project_root, "results", "tables", "11_TO_GCN", "input")
togcn_dir <- file.path(
    project_root, "results", "tables", "11_TO_GCN", "original_TO_GCN"
)
table_dir <- file.path(project_root, "results", "tables", "11_TO_GCN", "postprocess")
figure_dir <- file.path(project_root, "results", "figures", "11_TO_GCN")
session_dir <- file.path(project_root, "results", "sessionInfo")
reference_dir <- file.path(project_root, "data", "reference", "TO_GCN")

dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(session_dir, recursive = TRUE, showWarnings = FALSE)

domain_levels <- c("SAM", "P1_P2", "P3", "P4", "P5")
positive_pcc_cutoff <- 0.95
maximum_edges_to_plot <- 50000L

tf_mean_file <- file.path(input_dir, "TF_genes_mean_UMI_per_spot.csv")
tf_zscore_file <- file.path(input_dir, "TF_genes_gene_wise_zscore.csv")
all_gene_file <- file.path(input_dir, "all_genes_mean_UMI_per_spot.csv.gz")
annotation_file <- file.path(input_dir, "TF_annotation_used.csv")
node_level_file <- file.path(togcn_dir, "Node_level.csv")
node_relation_file <- file.path(togcn_dir, "Node_relation.csv")
tf_gene_edge_file <- file.path(togcn_dir, "C1+.csv")
reference_level_file <- file.path(input_dir, "recovered_reference_TF_levels.csv")
run_parameter_file <- file.path(togcn_dir, "run_parameters.tsv")

if (file.exists(run_parameter_file)) {
    run_parameters <- read.delim(
        run_parameter_file,
        stringsAsFactors = FALSE,
        check.names = FALSE
    )
    cutoff_row <- which(run_parameters$parameter == "positive_PCC_cutoff")
    if (length(cutoff_row) == 1L) {
        recorded_cutoff <- suppressWarnings(as.numeric(run_parameters$value[cutoff_row]))
        if (is.finite(recorded_cutoff)) positive_pcc_cutoff <- recorded_cutoff
    }
}

required_files <- c(tf_mean_file, tf_zscore_file, all_gene_file, annotation_file)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files) > 0L) {
    stop("Missing preprocessing output(s): ", paste(missing_files, collapse = ", "))
}

tf_mean <- read_expression_matrix(tf_mean_file)
tf_zscore <- read_expression_matrix(tf_zscore_file)
all_gene_mean <- read_expression_matrix(all_gene_file)
tf_annotation <- read.csv(
    annotation_file,
    stringsAsFactors = FALSE,
    check.names = FALSE
)

tf_mean <- tf_mean[, domain_levels, drop = FALSE]
tf_zscore <- tf_zscore[, domain_levels, drop = FALSE]

# Prefer the freshly generated TO-GCN levels. The recovered workbook levels are
# an explicit demonstration fallback, not a substitute for rerunning TO-GCN on
# a new dataset.
level_source <- "original TO-GCN Node_level.csv"
if (file.exists(node_level_file)) {
    levels_raw <- read.csv(
        node_level_file,
        stringsAsFactors = FALSE,
        check.names = FALSE
    )
    names(levels_raw) <- trimws(names(levels_raw))
    gene_column <- detect_column(names(levels_raw), c("tf_gene_id", "gene_id", "gene"))
    level_column <- detect_column(names(levels_raw), c("level_in_gcn", "level"))
    if (is.na(gene_column) || is.na(level_column)) {
        stop("Could not interpret columns in Node_level.csv.")
    }
    tf_levels <- data.frame(
        gene_id = trimws(as.character(levels_raw[[gene_column]])),
        level = suppressWarnings(as.integer(levels_raw[[level_column]])),
        stringsAsFactors = FALSE
    )
} else if (file.exists(reference_level_file)) {
    warning(
        "Node_level.csv is absent. Using recovered reference levels for the ",
        "demonstration figure. Rerun TO-GCN for a new dataset."
    )
    levels_raw <- read.csv(reference_level_file, stringsAsFactors = FALSE)
    tf_levels <- data.frame(
        gene_id = trimws(as.character(levels_raw$gene_id)),
        level = suppressWarnings(as.integer(levels_raw$recovered_level)),
        stringsAsFactors = FALSE
    )
    level_source <- "recovered reference workbook (demonstration fallback)"
} else {
    stop(
        "Neither Node_level.csv nor recovered_reference_TF_levels.csv exists. ",
        "Run scripts/shell/11_TO_GCN/02_run_TO_GCN_single.sh first."
    )
}

tf_levels <- tf_levels |>
    filter(!is.na(gene_id), nzchar(gene_id), !is.na(level), level > 0L) |>
    distinct(gene_id, .keep_all = TRUE) |>
    filter(gene_id %in% rownames(tf_zscore))

if (nrow(tf_levels) == 0L) stop("No assigned TF levels overlap the TF matrix.")

tf_level_annotation <- tf_levels |>
    left_join(tf_annotation, by = "gene_id") |>
    arrange(level, gene_id)

write.csv(
    tf_level_annotation,
    file.path(table_dir, "TF_level_assignments.csv"),
    row.names = FALSE,
    quote = FALSE
)

level_counts <- tf_level_annotation |>
    count(level, name = "n_TFs") |>
    arrange(level)
write.csv(
    level_counts,
    file.path(table_dir, "TF_counts_by_TO_GCN_level.csv"),
    row.names = FALSE,
    quote = FALSE
)

if (file.exists(reference_level_file) && file.exists(node_level_file)) {
    reference_levels <- read.csv(reference_level_file, stringsAsFactors = FALSE) |>
        transmute(
            gene_id = trimws(as.character(gene_id)),
            recovered_level = as.integer(recovered_level)
        )
    comparison <- tf_levels |>
        inner_join(reference_levels, by = "gene_id") |>
        mutate(level_identical = level == recovered_level)
    write.csv(
        comparison,
        file.path(table_dir, "new_vs_recovered_TF_level_comparison.csv"),
        row.names = FALSE,
        quote = FALSE
    )
    write.csv(
        data.frame(
            n_compared = nrow(comparison),
            n_identical = sum(comparison$level_identical),
            fraction_identical = mean(comparison$level_identical)
        ),
        file.path(table_dir, "new_vs_recovered_TF_level_summary.csv"),
        row.names = FALSE,
        quote = FALSE
    )
}

# -----------------------------
# Figure D: ordered TF heatmap
# -----------------------------

assigned_gene_ids <- intersect(tf_levels$gene_id, rownames(tf_zscore))
heatmap_matrix <- tf_zscore[assigned_gene_ids, domain_levels, drop = FALSE]
heatmap_levels <- tf_levels$level[match(rownames(heatmap_matrix), tf_levels$gene_id)]

# Order genes within each BFS level by their developmental center of expression.
# This affects only the display order, never the assigned TO-GCN level.
time_score <- as.vector(heatmap_matrix %*% seq_along(domain_levels))
display_order <- order(heatmap_levels, time_score, rownames(heatmap_matrix))
heatmap_matrix <- heatmap_matrix[display_order, , drop = FALSE]
heatmap_levels <- heatmap_levels[display_order]

gene_x <- setNames(seq_len(nrow(heatmap_matrix)), rownames(heatmap_matrix))
heatmap_long <- data.frame(
    gene_id = rep(rownames(heatmap_matrix), times = length(domain_levels)),
    domain = rep(domain_levels, each = nrow(heatmap_matrix)),
    z_score = as.vector(heatmap_matrix),
    stringsAsFactors = FALSE
)
heatmap_long$x <- unname(gene_x[heatmap_long$gene_id])
heatmap_long$domain <- factor(heatmap_long$domain, levels = rev(domain_levels))

level_layout <- data.frame(
    level = sort(unique(heatmap_levels)),
    start = vapply(sort(unique(heatmap_levels)), function(x) min(which(heatmap_levels == x)), numeric(1)),
    end = vapply(sort(unique(heatmap_levels)), function(x) max(which(heatmap_levels == x)), numeric(1))
) |>
    mutate(center = (start + end) / 2)

published_marker_names <- c(
    "KN1", "ARF5", "EREB158", "MYB15", "NS1", "ARF8", "ARF28",
    "MYB23", "GRAS25", "MYB85", "GRF15", "BHLH8", "EREB66",
    "GRAS54", "MYB11"
)
marker_annotations <- tf_level_annotation |>
    mutate(gene_name_upper = toupper(gene_name)) |>
    filter(gene_name_upper %in% published_marker_names, gene_id %in% names(gene_x)) |>
    mutate(x = unname(gene_x[gene_id])) |>
    arrange(match(gene_name_upper, published_marker_names)) |>
    distinct(gene_name_upper, .keep_all = TRUE)

panel_d <- ggplot(heatmap_long, aes(x = x, y = domain, fill = z_score)) +
    geom_raster() +
    geom_vline(
        data = level_layout[-nrow(level_layout), , drop = FALSE],
        aes(xintercept = end + 0.5),
        inherit.aes = FALSE,
        color = "white",
        linewidth = 0.25
    ) +
    geom_text(
        data = marker_annotations,
        aes(x = x, y = "SAM", label = gene_name_upper),
        inherit.aes = FALSE,
        angle = 90,
        vjust = -1.0,
        hjust = 0,
        size = 2.4,
        fontface = "italic"
    ) +
    scale_x_continuous(
        breaks = level_layout$center,
        labels = paste0("L", level_layout$level),
        expand = expansion(mult = c(0, 0))
    ) +
    scale_fill_gradient2(
        low = "#172A88",
        mid = "white",
        high = "#B2182B",
        midpoint = 0,
        limits = c(-1.5, 1.5),
        oob = clamp_values,
        name = "Z-score"
    ) +
    coord_cartesian(clip = "off") +
    labs(x = NULL, y = NULL, title = "D") +
    theme_classic(base_size = 11) +
    theme(
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_text(size = 9),
        axis.text.y = element_text(face = "bold"),
        plot.title = element_text(face = "bold", size = 16),
        plot.margin = margin(t = 35, r = 5, b = 5, l = 5),
        legend.position = "right"
    )

# -----------------------------
# Figure E: level network and mean expression profiles
# -----------------------------

if (file.exists(node_relation_file)) {
    relation_raw <- read.csv(
        node_relation_file,
        stringsAsFactors = FALSE,
        check.names = FALSE
    )
    names(relation_raw) <- trimws(names(relation_raw))
    node1_column <- detect_column(names(relation_raw), c("node_1_id", "node1_id"))
    node2_column <- detect_column(names(relation_raw), c("node_2_id", "node2_id"))
    pcc_column <- detect_column(names(relation_raw), c("pcc_value", "pcc"))
    if (is.na(node1_column) || is.na(node2_column)) {
        stop("Could not interpret columns in Node_relation.csv.")
    }
    network_edges <- data.frame(
        from = trimws(as.character(relation_raw[[node1_column]])),
        to = trimws(as.character(relation_raw[[node2_column]])),
        PCC = if (is.na(pcc_column)) NA_real_ else
            suppressWarnings(as.numeric(relation_raw[[pcc_column]])),
        stringsAsFactors = FALSE
    )
} else {
    warning("Node_relation.csv is absent; recalculating TF-TF PCC edges for plotting.")
    common_tf <- intersect(tf_levels$gene_id, rownames(tf_mean))
    correlation_matrix <- cor(t(tf_mean[common_tf, , drop = FALSE]))
    edge_indices <- which(
        upper.tri(correlation_matrix) & correlation_matrix >= positive_pcc_cutoff,
        arr.ind = TRUE
    )
    network_edges <- data.frame(
        from = rownames(correlation_matrix)[edge_indices[, 1L]],
        to = colnames(correlation_matrix)[edge_indices[, 2L]],
        PCC = correlation_matrix[edge_indices],
        stringsAsFactors = FALSE
    )
}

network_edges <- network_edges |>
    filter(from %in% tf_levels$gene_id, to %in% tf_levels$gene_id) |>
    distinct(from, to, .keep_all = TRUE)

node_positions <- bind_rows(lapply(sort(unique(tf_levels$level)), function(current_level) {
    genes <- sort(tf_levels$gene_id[tf_levels$level == current_level])
    n_nodes <- length(genes)
    maximum_nodes <- max(level_counts$n_TFs)
    radius <- 0.22 + 0.48 * sqrt(n_nodes / maximum_nodes)
    angle <- if (n_nodes == 1L) 0 else seq(0, 2 * pi, length.out = n_nodes + 1L)[-1L]
    data.frame(
        gene_id = genes,
        level = current_level,
        center_x = current_level * 1.45,
        radius = radius,
        x = current_level * 1.45 + radius * cos(angle),
        y = radius * sin(angle),
        stringsAsFactors = FALSE
    )
}))

if (nrow(network_edges) > maximum_edges_to_plot) {
    network_edges <- network_edges[
        sample(seq_len(nrow(network_edges)), maximum_edges_to_plot),
        ,
        drop = FALSE
    ]
}

edge_plot_data <- network_edges |>
    left_join(
        node_positions |> select(from = gene_id, x, y),
        by = "from"
    ) |>
    rename(x_from = x, y_from = y) |>
    left_join(
        node_positions |> select(to = gene_id, x, y),
        by = "to"
    ) |>
    rename(x_to = x, y_to = y) |>
    filter(complete.cases(x_from, y_from, x_to, y_to))

circle_angles <- seq(0, 2 * pi, length.out = 201L)
circle_data <- bind_rows(lapply(split(node_positions, node_positions$level), function(x) {
    data.frame(
        level = x$level[1L],
        x = x$center_x[1L] + x$radius[1L] * cos(circle_angles),
        y = x$radius[1L] * sin(circle_angles)
    )
}))

network_labels <- node_positions |>
    group_by(level) |>
    summarise(
        x = first(center_x),
        radius = first(radius),
        n_TFs = n(),
        .groups = "drop"
    )

network_plot <- ggplot() +
    geom_segment(
        data = edge_plot_data,
        aes(x = x_from, y = y_from, xend = x_to, yend = y_to),
        color = "#91A9C9",
        alpha = 0.035,
        linewidth = 0.2
    ) +
    geom_path(
        data = circle_data,
        aes(x = x, y = y, group = level),
        color = "#173C88",
        linewidth = 0.55
    ) +
    geom_point(
        data = node_positions,
        aes(x = x, y = y),
        color = "#173C88",
        size = 0.22
    ) +
    geom_text(
        data = network_labels,
        aes(x = x, y = 0, label = n_TFs),
        size = 3.1
    ) +
    geom_text(
        data = network_labels,
        aes(x = x, y = -radius - 0.18, label = paste0("L", level)),
        fontface = "bold",
        size = 3.2
    ) +
    coord_equal(clip = "off") +
    labs(title = "E", subtitle = paste0("TO-GCN TF network; level source: ", level_source)) +
    theme_void(base_size = 11) +
    theme(
        plot.title = element_text(face = "bold", size = 16),
        plot.subtitle = element_text(size = 8),
        plot.margin = margin(5, 10, 15, 10)
    )

level_profiles <- tf_level_annotation |>
    filter(gene_id %in% rownames(tf_zscore)) |>
    select(gene_id, level) |>
    left_join(
        data.frame(
            gene_id = rownames(tf_zscore),
            tf_zscore,
            check.names = FALSE
        ),
        by = "gene_id"
    ) |>
    pivot_longer(
        cols = all_of(domain_levels),
        names_to = "domain",
        values_to = "z_score"
    ) |>
    group_by(level, domain) |>
    summarise(mean_z_score = mean(z_score), .groups = "drop") |>
    mutate(
        domain = factor(domain, levels = domain_levels),
        level_label = factor(
            paste0("L", level),
            levels = paste0("L", sort(unique(level)))
        ),
        sign = if_else(mean_z_score >= 0, "Positive", "Negative")
    )

write.csv(
    level_profiles,
    file.path(table_dir, "mean_zscore_profile_by_TO_GCN_level.csv"),
    row.names = FALSE,
    quote = FALSE
)

profile_plot <- ggplot(
    level_profiles,
    aes(x = domain, y = mean_z_score, fill = sign)
) +
    geom_hline(yintercept = 0, linewidth = 0.25) +
    geom_col(width = 0.72) +
    facet_wrap(~level_label, ncol = 7) +
    scale_fill_manual(values = c(Negative = "#2456C4", Positive = "#E53B25")) +
    labs(x = NULL, y = "Mean TF z-score", fill = NULL) +
    theme_classic(base_size = 8) +
    theme(
        strip.background = element_blank(),
        strip.text = element_text(face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 6),
        legend.position = "bottom"
    )

panel_e <- network_plot / profile_plot + plot_layout(heights = c(1.15, 1.6))
combined_figure <- panel_d / panel_e + plot_layout(heights = c(1.05, 1.65))

ggsave(
    file.path(figure_dir, "Figure_TO_GCN_D_TF_heatmap.png"),
    panel_d,
    width = 14,
    height = 4.8,
    dpi = 400,
    bg = "white"
)
ggsave(
    file.path(figure_dir, "Figure_TO_GCN_E_level_network_and_profiles.png"),
    panel_e,
    width = 14,
    height = 7.0,
    dpi = 400,
    bg = "white"
)
ggsave(
    file.path(figure_dir, "Figure_TO_GCN_D_E_composite.png"),
    combined_figure,
    width = 14,
    height = 11.5,
    dpi = 400,
    bg = "white"
)
ggsave(
    file.path(figure_dir, "Figure_TO_GCN_D_E_composite.pdf"),
    combined_figure,
    width = 14,
    height = 11.5,
    device = "pdf",
    bg = "white"
)

# -----------------------------
# Level-specific TF + coexpressed-gene sets and optional GO enrichment
# -----------------------------

background_genes <- rownames(all_gene_mean)
writeLines(
    background_genes,
    file.path(table_dir, "all_expressed_genes_AgriGO_background.txt")
)

if (file.exists(tf_gene_edge_file)) {
    tf_gene_edges_raw <- read.csv(
        tf_gene_edge_file,
        stringsAsFactors = FALSE,
        check.names = FALSE
    )
    names(tf_gene_edges_raw) <- trimws(names(tf_gene_edges_raw))
    tf_column <- detect_column(names(tf_gene_edges_raw), c("tf_gene_id"))
    target_column <- detect_column(names(tf_gene_edges_raw), c("gene_id"))
    pcc_column <- detect_column(names(tf_gene_edges_raw), c("pcc_under_c1", "pcc"))
    if (is.na(tf_column) || is.na(target_column)) {
        stop("Could not interpret columns in C1+.csv.")
    }
    tf_gene_edges <- data.frame(
        tf_gene_id = trimws(as.character(tf_gene_edges_raw[[tf_column]])),
        target_gene_id = trimws(as.character(tf_gene_edges_raw[[target_column]])),
        PCC = if (is.na(pcc_column)) NA_real_ else
            suppressWarnings(as.numeric(tf_gene_edges_raw[[pcc_column]])),
        stringsAsFactors = FALSE
    ) |>
        filter(tf_gene_id %in% tf_levels$gene_id) |>
        left_join(tf_levels, by = c("tf_gene_id" = "gene_id")) |>
        filter(!is.na(level))
} else {
    warning(
        "C1+.csv is absent. AgriGO lists will contain TF genes only. ",
        "Run GCN_single to include coexpressed non-TF targets."
    )
    tf_gene_edges <- data.frame(
        tf_gene_id = character(), target_gene_id = character(),
        PCC = numeric(), level = integer()
    )
}

level_gene_sets <- lapply(sort(unique(tf_levels$level)), function(current_level) {
    unique(c(
        tf_levels$gene_id[tf_levels$level == current_level],
        tf_gene_edges$target_gene_id[tf_gene_edges$level == current_level]
    ))
})
names(level_gene_sets) <- paste0("L", sort(unique(tf_levels$level)))
level_gene_sets <- lapply(level_gene_sets, intersect, y = background_genes)

level_set_summary <- data.frame(
    level = names(level_gene_sets),
    n_genes = lengths(level_gene_sets),
    stringsAsFactors = FALSE
)
write.csv(
    level_set_summary,
    file.path(table_dir, "TO_GCN_level_gene_set_sizes.csv"),
    row.names = FALSE,
    quote = FALSE
)

for (level_name in names(level_gene_sets)) {
    writeLines(
        level_gene_sets[[level_name]],
        file.path(table_dir, paste0(level_name, "_TF_and_coexpressed_genes_AgriGO.txt"))
    )
}

gene_to_go_file <- file.path(reference_dir, "maize_gene_to_GO.csv")
if (file.exists(gene_to_go_file)) {
    gene_to_go_raw <- read.csv(
        gene_to_go_file,
        stringsAsFactors = FALSE,
        check.names = FALSE
    )
    gene_column <- detect_column(names(gene_to_go_raw), c("gene_id", "gene"))
    go_column <- detect_column(names(gene_to_go_raw), c("go_id", "go"))
    if (is.na(gene_column) || is.na(go_column)) {
        stop("maize_gene_to_GO.csv requires gene_id and GO_ID columns.")
    }
    gene_to_go <- data.frame(
        gene_id = trimws(as.character(gene_to_go_raw[[gene_column]])),
        GO_ID = as.character(gene_to_go_raw[[go_column]]),
        stringsAsFactors = FALSE
    ) |>
        separate_rows(GO_ID, sep = "[;,|[:space:]]+") |>
        mutate(GO_ID = toupper(trimws(GO_ID))) |>
        filter(gene_id %in% background_genes, grepl("^GO:[0-9]+$", GO_ID)) |>
        distinct()

    annotated_background <- intersect(background_genes, unique(gene_to_go$gene_id))
    background_size <- length(annotated_background)
    go_sizes <- table(gene_to_go$GO_ID)

    enrichment <- bind_rows(lapply(names(level_gene_sets), function(level_name) {
        query <- intersect(level_gene_sets[[level_name]], annotated_background)
        if (length(query) == 0L) return(NULL)
        overlap <- table(gene_to_go$GO_ID[gene_to_go$gene_id %in% query])
        go_ids <- names(overlap)
        term_sizes <- as.integer(go_sizes[go_ids])
        overlap_sizes <- as.integer(overlap)
        data.frame(
            level = level_name,
            GO_ID = go_ids,
            overlap = overlap_sizes,
            query_size = length(query),
            term_size = term_sizes,
            background_size = background_size,
            p_value = phyper(
                overlap_sizes - 1L,
                term_sizes,
                background_size - term_sizes,
                length(query),
                lower.tail = FALSE
            ),
            stringsAsFactors = FALSE
        ) |>
            mutate(FDR = p.adjust(p_value, method = "BH"))
    }))

    write.csv(
        enrichment,
        file.path(table_dir, "TO_GCN_level_GO_Fisher_BH_enrichment.csv"),
        row.names = FALSE,
        quote = FALSE
    )
}

write.csv(
    data.frame(
        parameter = c(
            "level_source", "ordered_domains", "positive_PCC_cutoff",
            "maximum_edges_plotted", "GO_test", "multiple_testing"
        ),
        value = c(
            level_source, paste(domain_levels, collapse = ";"),
            positive_pcc_cutoff, maximum_edges_to_plot,
            "one-sided Fisher exact / hypergeometric", "Benjamini-Hochberg FDR"
        )
    ),
    file.path(table_dir, "TO_GCN_postprocess_parameters.csv"),
    row.names = FALSE,
    quote = TRUE
)

writeLines(
    capture.output(sessionInfo()),
    file.path(session_dir, "11_TO_GCN_postprocess_sessionInfo.txt")
)

message("TO-GCN postprocessing completed.")
message("Level source: ", level_source)
message("Assigned TFs: ", nrow(tf_levels))
message("Levels: ", paste(sort(unique(tf_levels$level)), collapse = ", "))
message("Figures: ", figure_dir)
