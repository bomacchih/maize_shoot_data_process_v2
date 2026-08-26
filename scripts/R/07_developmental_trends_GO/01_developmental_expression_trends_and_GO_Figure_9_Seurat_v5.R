# Developmental expression-trend clustering and GO enrichment (Figure 9)
#
# Project input:
#   data/processed/maize_shoot_14samples_SCT_harmony_seurat_v5.rds
#
# Optional reference inputs:
#   data/reference/developmental_trends/sub_gene.csv
#   data/reference/developmental_trends/zea_go2.csv
#   data/reference/developmental_trends/maize_id_name.csv
#   data/reference/developmental_trends/go_term_descriptions.csv
#
# Outputs:
#   results/tables/07_developmental_trends_GO/
#   results/figures/07_developmental_trends_GO/
#   results/sessionInfo/07_developmental_trends_GO_sessionInfo.txt

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
})

# -----------------------------
# 1. Paths and analysis settings
# -----------------------------

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

project_root <- find_project_root()

input_rds <- file.path(
  project_root,
  "data", "processed",
  "maize_shoot_14samples_SCT_harmony_seurat_v5.rds"
)

reference_dir <- file.path(
  project_root, "data", "reference", "developmental_trends"
)
table_dir <- file.path(
  project_root, "results", "tables", "07_developmental_trends_GO"
)
figure_dir <- file.path(
  project_root, "results", "figures", "07_developmental_trends_GO"
)
session_dir <- file.path(project_root, "results", "sessionInfo")

dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(session_dir, recursive = TRUE, showWarnings = FALSE)

domain_levels <- c(
  "SAM", "P1_P2", "P3", "P4", "P5", "co_v", "coleoptile"
)

replicate_column <- "sample_id"
domain_column <- "domains"
number_of_clusters <- 7L
representatives_per_cluster <- 3L
go_genes_per_cluster <- 1000L
go_panel_clusters <- c("C1", "C6", "C7")
go_pvalue_cutoff <- 1e-6
go_terms_to_display <- 29L

# The seven-cluster cut was selected in the original study as the most
# informative summary of developmental expression patterns. The top-1,000-gene
# cutoff for GO analysis was an explicitly defined analytical choice.

# -----------------------------
# 2. Helper functions
# -----------------------------

get_assay_data_v5 <- function(object, assay = "SCT", layer = "data") {
  if (!assay %in% Assays(object)) {
    stop("Assay '", assay, "' is absent. Available assays: ",
         paste(Assays(object), collapse = ", "))
  }

  # Match layer names locally. Passing the anchored expression to Layers(search=)
  # is unreliable for SCTAssay objects in some SeuratObject releases because it
  # may be interpreted as a literal layer name.
  all_layer_names <- Layers(object[[assay]])
  layer_names <- grep(
    paste0("^", layer, "($|\\.)"),
    all_layer_names,
    value = TRUE
  )
  if (length(layer_names) == 0L) {
    stop("No '", layer, "' layer was found in the ", assay, " assay.")
  }

  matrices <- lapply(layer_names, function(layer_name) {
    LayerData(object, assay = assay, layer = layer_name)
  })

  if (length(matrices) == 1L) {
    return(matrices[[1L]])
  }

  common_features <- Reduce(intersect, lapply(matrices, rownames))
  if (length(common_features) == 0L) {
    stop("The SCT data layers do not share any features.")
  }

  matrices <- lapply(
    matrices,
    function(x) x[common_features, , drop = FALSE]
  )
  combined <- do.call(cbind, matrices)

  if (anyDuplicated(colnames(combined))) {
    stop("Duplicated spot names were found while combining SCT data layers.")
  }
  combined
}

row_scale <- function(x) {
  x <- as.matrix(x)
  gene_mean <- rowMeans(x)
  centered <- sweep(x, 1L, gene_mean, FUN = "-")
  gene_sd <- sqrt(rowSums(centered^2) / (ncol(x) - 1L))
  keep <- is.finite(gene_sd) & gene_sd > 0
  scaled <- sweep(centered[keep, , drop = FALSE], 1L, gene_sd[keep], FUN = "/")
  list(matrix = scaled, retained = rownames(x)[keep])
}

detect_column <- function(column_names, patterns, fallback = NA_character_) {
  normalized <- tolower(gsub("[^a-z0-9]+", "_", column_names))
  for (pattern in patterns) {
    hit <- grep(pattern, normalized, perl = TRUE)
    if (length(hit) > 0L) return(column_names[hit[1L]])
  }
  fallback
}

read_gene_list <- function(path) {
  x <- read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  fallback_col <- names(x)[1L]
  if (ncol(x) > 1L &&
      (is.na(fallback_col) || fallback_col == "" || grepl("^X(\\.[0-9]+)?$", fallback_col))) {
    fallback_col <- names(x)[2L]
  }
  gene_col <- detect_column(
    names(x),
    c("^gene$", "gene_id", "locus", "feature", "^x$"),
    fallback_col
  )
  unique(na.omit(trimws(as.character(x[[gene_col]]))))
}

read_gene_symbols <- function(path) {
  x <- read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  gene_col <- detect_column(names(x), c("gene_id", "^gene$", "locus", "feature"), names(x)[1L])
  symbol_col <- detect_column(names(x), c("symbol", "gene_name", "short_name", "name"), NA_character_)
  if (is.na(symbol_col)) return(setNames(character(), character()))
  x <- x[!is.na(x[[gene_col]]) & !is.na(x[[symbol_col]]), , drop = FALSE]
  setNames(as.character(x[[symbol_col]]), as.character(x[[gene_col]]))
}

read_gene_to_go <- function(path) {
  x <- read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  gene_col <- detect_column(names(x), c("gene_id", "^gene$", "locus", "feature"), names(x)[1L])
  go_col <- detect_column(names(x), c("go_id", "^go$", "gene_ontology"), names(x)[min(2L, ncol(x))])

  out <- data.frame(
    gene = trimws(as.character(x[[gene_col]])),
    go_entry = as.character(x[[go_col]]),
    stringsAsFactors = FALSE
  ) |>
    filter(!is.na(gene), gene != "", !is.na(go_entry), go_entry != "") |>
    separate_rows(go_entry, sep = "[;,|[:space:]]+") |>
    transmute(gene, GO_ID = toupper(trimws(go_entry))) |>
    filter(grepl("^GO:[0-9]+$", GO_ID)) |>
    distinct()

  if (nrow(out) == 0L) {
    stop("No valid GO identifiers were parsed from ", path)
  }
  out
}

hypergeometric_go <- function(query_by_cluster, gene_to_go, background_genes) {
  gene_to_go <- gene_to_go |>
    filter(gene %in% background_genes) |>
    distinct(gene, GO_ID)

  annotated_background <- intersect(background_genes, unique(gene_to_go$gene))
  total_background <- length(annotated_background)
  go_sizes <- table(gene_to_go$GO_ID)

  results <- lapply(names(query_by_cluster), function(cluster_name) {
    query <- intersect(unique(query_by_cluster[[cluster_name]]), annotated_background)
    if (length(query) == 0L) return(NULL)

    query_map <- gene_to_go |> filter(gene %in% query)
    overlap <- table(query_map$GO_ID)
    go_ids <- names(overlap)
    term_size <- as.integer(go_sizes[go_ids])
    overlap_size <- as.integer(overlap)

    data.frame(
      cluster = cluster_name,
      GO_ID = go_ids,
      overlap = overlap_size,
      query_size = length(query),
      term_size = term_size,
      background_size = total_background,
      p_value = phyper(
        overlap_size - 1L,
        term_size,
        total_background - term_size,
        length(query),
        lower.tail = FALSE
      ),
      stringsAsFactors = FALSE
    ) |>
      mutate(p_adj = p.adjust(p_value, method = "BH")) |>
      arrange(p_value)
  })

  bind_rows(results)
}

add_go_descriptions <- function(results, reference_directory) {
  description_file <- file.path(reference_directory, "go_term_descriptions.csv")
  term_map <- NULL

  if (file.exists(description_file)) {
    x <- read.csv(description_file, check.names = FALSE, stringsAsFactors = FALSE)
    go_col <- detect_column(names(x), c("go_id", "^go$"), names(x)[1L])
    term_col <- detect_column(
      names(x), c("description", "^term$", "go_term", "name"),
      names(x)[min(2L, ncol(x))]
    )
    term_map <- data.frame(
      GO_ID = as.character(x[[go_col]]),
      GO_term = as.character(x[[term_col]]),
      stringsAsFactors = FALSE
    ) |>
      distinct(GO_ID, .keep_all = TRUE)
  } else if (requireNamespace("GO.db", quietly = TRUE) &&
             requireNamespace("AnnotationDbi", quietly = TRUE)) {
    term_map <- AnnotationDbi::select(
      GO.db::GO.db,
      keys = unique(results$GO_ID),
      keytype = "GOID",
      columns = "TERM"
    ) |>
      transmute(GO_ID = GOID, GO_term = TERM) |>
      distinct(GO_ID, .keep_all = TRUE)
  }

  if (is.null(term_map)) {
    results$GO_term <- results$GO_ID
  } else {
    results <- results |>
      left_join(term_map, by = "GO_ID") |>
      mutate(GO_term = if_else(is.na(GO_term) | GO_term == "", GO_ID, GO_term))
  }
  results
}

# -----------------------------
# 3. Load data and aggregate by biological replicate and domain
# -----------------------------

if (!file.exists(input_rds)) stop("Input RDS not found: ", input_rds)
maize <- readRDS(input_rds)

required_metadata <- c(replicate_column, domain_column)
missing_metadata <- setdiff(required_metadata, colnames(maize[[]]))
if (length(missing_metadata) > 0L) {
  stop("Missing metadata column(s): ", paste(missing_metadata, collapse = ", "))
}

metadata <- maize[[]]
metadata$spot <- rownames(metadata)
metadata$replicate <- as.character(metadata[[replicate_column]])
metadata$domain <- as.character(metadata[[domain_column]])
metadata <- metadata |>
  filter(!is.na(replicate), replicate != "", domain %in% domain_levels)

if (nrow(metadata) == 0L) stop("No spots remained in the seven target domains.")
metadata$domain <- factor(metadata$domain, levels = domain_levels)

sct_data <- get_assay_data_v5(maize, assay = "SCT", layer = "data")
missing_spots <- setdiff(metadata$spot, colnames(sct_data))
if (length(missing_spots) > 0L) {
  stop(length(missing_spots), " metadata spot(s) are absent from the SCT data layer.")
}
sct_data <- sct_data[, metadata$spot, drop = FALSE]

group_id <- interaction(
  metadata$replicate,
  metadata$domain,
  sep = "__",
  drop = TRUE,
  lex.order = TRUE
)
group_levels <- levels(group_id)
membership <- sparseMatrix(
  i = seq_along(group_id),
  j = as.integer(group_id),
  x = 1,
  dims = c(length(group_id), length(group_levels)),
  dimnames = list(metadata$spot, group_levels)
)

group_spot_counts <- Matrix::colSums(membership)
replicate_domain_expression <- sct_data %*% membership
replicate_domain_expression <- sweep(
  replicate_domain_expression,
  2L,
  group_spot_counts,
  FUN = "/"
)

group_metadata <- data.frame(
  group = group_levels,
  replicate = sub("__.*$", "", group_levels),
  domain = sub("^.*__", "", group_levels),
  n_spots = as.integer(group_spot_counts),
  stringsAsFactors = FALSE
)
group_metadata$domain <- factor(group_metadata$domain, levels = domain_levels)

write.csv(
  group_metadata,
  file.path(table_dir, "replicate_domain_profiles.csv"),
  row.names = FALSE
)

# Retain genes detected in every target domain. Detection is assessed from the
# replicate-averaged SCT expression and does not treat individual spots as
# independent biological replicates.
domain_unscaled <- sapply(domain_levels, function(domain_name) {
  cols <- which(group_metadata$domain == domain_name)
  if (length(cols) == 0L) stop("No replicate profile was available for ", domain_name)
  Matrix::rowMeans(replicate_domain_expression[, cols, drop = FALSE])
})

detected_in_all_domains <- rowSums(domain_unscaled != 0) == length(domain_levels)
eligible_genes <- rownames(domain_unscaled)[detected_in_all_domains]

gene_list_file <- file.path(reference_dir, "sub_gene.csv")
if (file.exists(gene_list_file)) {
  recovered_gene_list <- read_gene_list(gene_list_file)
  clustering_genes <- intersect(eligible_genes, recovered_gene_list)
  message("Using ", length(clustering_genes),
          " eligible genes shared with the recovered gene list.")
} else {
  clustering_genes <- eligible_genes
  message("Recovered sub_gene.csv was not found; using all ",
          length(clustering_genes), " eligible genes.")
}

if (length(clustering_genes) < number_of_clusters) {
  stop("Too few eligible genes for seven-cluster hierarchical clustering.")
}

# Scale within each gene across replicate-by-domain profiles, then average the
# scaled profiles across biological replicates for descriptive domain profiles.
scaled_result <- row_scale(replicate_domain_expression[clustering_genes, , drop = FALSE])
replicate_domain_scaled <- scaled_result$matrix
clustering_genes <- scaled_result$retained

domain_mean_scaled <- sapply(domain_levels, function(domain_name) {
  cols <- which(group_metadata$domain == domain_name)
  rowMeans(replicate_domain_scaled[, cols, drop = FALSE])
})
colnames(domain_mean_scaled) <- domain_levels

write.csv(
  data.frame(gene = rownames(domain_mean_scaled), domain_mean_scaled, check.names = FALSE),
  file.path(table_dir, "domain_mean_scaled_expression.csv"),
  row.names = FALSE
)

# -----------------------------
# 4. Hierarchical clustering into seven expression trends
# -----------------------------

gene_distance <- dist(domain_mean_scaled, method = "euclidean")
gene_tree <- hclust(gene_distance, method = "complete")
raw_cluster <- cutree(gene_tree, k = number_of_clusters)

# cutree labels have no biological meaning. Relabel clusters from early to late
# expression using the developmental center of each cluster median. This changes
# labels only; it does not change cluster membership.
raw_medians <- lapply(sort(unique(raw_cluster)), function(cluster_id) {
  apply(domain_mean_scaled[raw_cluster == cluster_id, , drop = FALSE], 2L, median)
})
names(raw_medians) <- sort(unique(raw_cluster))

developmental_center <- vapply(raw_medians, function(profile) {
  weight <- profile - min(profile) + 1e-8
  sum(seq_along(profile) * weight) / sum(weight)
}, numeric(1))

ordered_raw_clusters <- names(sort(developmental_center))
cluster_map <- setNames(
  paste0("C", seq_along(ordered_raw_clusters)),
  ordered_raw_clusters
)
trend_cluster <- unname(cluster_map[as.character(raw_cluster)])
trend_cluster <- factor(trend_cluster, levels = paste0("C", 1:number_of_clusters))

cluster_assignments <- data.frame(
  gene = rownames(domain_mean_scaled),
  raw_cutree_cluster = as.integer(raw_cluster),
  trend_cluster = trend_cluster,
  stringsAsFactors = FALSE
)

write.csv(
  cluster_assignments,
  file.path(table_dir, "seven_expression_trend_cluster_assignments.csv"),
  row.names = FALSE
)
saveRDS(
  gene_tree,
  file.path(table_dir, "seven_expression_trend_hclust.rds")
)

# -----------------------------
# 5. Representative genes and Figure 9A
# -----------------------------

cluster_medians <- bind_rows(lapply(levels(trend_cluster), function(cluster_name) {
  genes <- cluster_assignments$gene[cluster_assignments$trend_cluster == cluster_name]
  profile <- apply(domain_mean_scaled[genes, , drop = FALSE], 2L, median)
  data.frame(
    trend_cluster = cluster_name,
    domain = factor(domain_levels, levels = domain_levels),
    median_scaled_expression = as.numeric(profile)
  )
}))

representative_genes <- bind_rows(lapply(levels(trend_cluster), function(cluster_name) {
  genes <- cluster_assignments$gene[cluster_assignments$trend_cluster == cluster_name]
  profiles <- domain_mean_scaled[genes, , drop = FALSE]
  median_profile <- apply(profiles, 2L, median)
  distances <- sqrt(rowSums(sweep(profiles, 2L, median_profile, FUN = "-")^2))
  selected <- names(sort(distances))[seq_len(min(representatives_per_cluster, length(distances)))]
  data.frame(
    trend_cluster = cluster_name,
    gene = selected,
    distance_to_cluster_median = unname(distances[selected]),
    stringsAsFactors = FALSE
  )
}))

symbol_file <- file.path(reference_dir, "maize_id_name.csv")
gene_symbols <- if (file.exists(symbol_file)) read_gene_symbols(symbol_file) else setNames(character(), character())
representative_genes$label <- unname(gene_symbols[representative_genes$gene])
representative_genes$label[
  is.na(representative_genes$label) | representative_genes$label == ""
] <- representative_genes$gene[
  is.na(representative_genes$label) | representative_genes$label == ""
]

write.csv(
  representative_genes,
  file.path(table_dir, "representative_genes_three_per_cluster.csv"),
  row.names = FALSE
)

trend_long <- data.frame(
  gene = rep(rownames(domain_mean_scaled), times = ncol(domain_mean_scaled)),
  domain = rep(domain_levels, each = nrow(domain_mean_scaled)),
  scaled_expression = as.vector(domain_mean_scaled),
  stringsAsFactors = FALSE
) |>
  left_join(cluster_assignments[, c("gene", "trend_cluster")], by = "gene") |>
  mutate(
    domain = factor(domain, levels = domain_levels),
    trend_cluster = factor(trend_cluster, levels = paste0("C", 1:number_of_clusters))
  )

label_positions <- trend_long |>
  group_by(trend_cluster) |>
  summarise(y_max = max(scaled_expression, na.rm = TRUE), .groups = "drop") |>
  left_join(
    representative_genes |>
      group_by(trend_cluster) |>
      arrange(distance_to_cluster_median, .by_group = TRUE) |>
      mutate(label_rank = row_number()) |>
      ungroup(),
    by = "trend_cluster"
  ) |>
  mutate(
    domain = factor("SAM", levels = domain_levels),
    y = y_max - (label_rank - 1L) * 0.28
  )

panel_a <- ggplot(trend_long, aes(domain, scaled_expression, group = gene)) +
  geom_line(color = "grey55", alpha = 0.18, linewidth = 0.18) +
  geom_line(
    data = cluster_medians,
    aes(y = median_scaled_expression, group = trend_cluster),
    color = "red", linewidth = 1.0
  ) +
  geom_text(
    data = label_positions,
    aes(x = domain, y = y, label = label),
    inherit.aes = FALSE,
    hjust = 0, fontface = "italic", size = 2.7
  ) +
  facet_wrap(~trend_cluster, nrow = 1, scales = "free_y") +
  labs(
    x = "Structural domains",
    y = "Mean scaled SCT expression"
  ) +
  theme_bw(base_size = 10) +
  theme(
    strip.background = element_rect(fill = "grey75"),
    strip.text = element_text(face = "bold", size = 11),
    axis.text.x = element_text(angle = 50, hjust = 1),
    panel.grid.minor = element_blank()
  )

ggsave(
  file.path(figure_dir, "Figure_9A_expression_trend_clusters.png"),
  panel_a, width = 13.5, height = 5.2, dpi = 300
)
ggsave(
  file.path(figure_dir, "Figure_9A_expression_trend_clusters.pdf"),
  panel_a, width = 13.5, height = 5.2
)

# -----------------------------
# 6. GO gene lists, enrichment, and Figure 9B
# -----------------------------

mean_expression <- Matrix::rowMeans(
  replicate_domain_expression[cluster_assignments$gene, , drop = FALSE]
)

go_gene_table <- cluster_assignments |>
  mutate(mean_SCT_expression = unname(mean_expression[gene])) |>
  group_by(trend_cluster) |>
  arrange(desc(mean_SCT_expression), gene, .by_group = TRUE) |>
  slice_head(n = go_genes_per_cluster) |>
  ungroup()

write.csv(
  go_gene_table,
  file.path(table_dir, "top_1000_genes_per_cluster_for_GO.csv"),
  row.names = FALSE
)

for (cluster_name in levels(trend_cluster)) {
  writeLines(
    go_gene_table$gene[go_gene_table$trend_cluster == cluster_name],
    file.path(table_dir, paste0(cluster_name, "_top1000_AgriGO_gene_list.txt"))
  )
}

go_map_file <- file.path(reference_dir, "zea_go2.csv")
if (file.exists(go_map_file)) {
  gene_to_go <- read_gene_to_go(go_map_file)
  query_by_cluster <- split(go_gene_table$gene, go_gene_table$trend_cluster)
  go_results <- hypergeometric_go(
    query_by_cluster = query_by_cluster,
    gene_to_go = gene_to_go,
    background_genes = eligible_genes
  ) |>
    add_go_descriptions(reference_directory = reference_dir) |>
    mutate(minus_log10_p = -log10(pmax(p_value, .Machine$double.xmin)))

  write.csv(
    go_results,
    file.path(table_dir, "GO_hypergeometric_enrichment_all_clusters.csv"),
    row.names = FALSE
  )

  selected_terms <- go_results |>
    filter(cluster %in% go_panel_clusters, p_value < go_pvalue_cutoff) |>
    group_by(GO_ID) |>
    summarise(best_p = min(p_value), .groups = "drop") |>
    arrange(best_p) |>
    slice_head(n = go_terms_to_display) |>
    pull(GO_ID)

  if (length(selected_terms) > 0L) {
    heatmap_data <- expand_grid(
      cluster = factor(go_panel_clusters, levels = rev(go_panel_clusters)),
      GO_ID = selected_terms
    ) |>
      left_join(
        go_results |>
          filter(cluster %in% go_panel_clusters, GO_ID %in% selected_terms) |>
          select(cluster, GO_ID, GO_term, minus_log10_p),
        by = c("cluster", "GO_ID")
      ) |>
      mutate(
        minus_log10_p = replace_na(minus_log10_p, 0),
        GO_term = if_else(is.na(GO_term) | GO_term == "", as.character(GO_ID), GO_term)
      )

    go_matrix <- heatmap_data |>
      select(cluster, GO_ID, minus_log10_p) |>
      pivot_wider(names_from = GO_ID, values_from = minus_log10_p) |>
      as.data.frame()
    rownames(go_matrix) <- as.character(go_matrix$cluster)
    go_matrix$cluster <- NULL

    if (ncol(go_matrix) > 1L) {
      go_order <- colnames(go_matrix)[hclust(dist(t(as.matrix(go_matrix))), method = "complete")$order]
    } else {
      go_order <- colnames(go_matrix)
    }

    label_map <- heatmap_data |>
      distinct(GO_ID, GO_term) |>
      filter(GO_ID %in% go_order)
    ordered_labels <- label_map$GO_term[match(go_order, label_map$GO_ID)]
    heatmap_data$GO_term <- factor(heatmap_data$GO_term, levels = ordered_labels)
    panel_b <- ggplot(heatmap_data, aes(GO_term, cluster, fill = minus_log10_p)) +
      geom_tile(color = "grey75", linewidth = 0.25) +
      scale_fill_gradient(low = "white", high = "red") +
      labs(
        x = "GO term",
        y = NULL,
        fill = expression(-log[10](italic(p))),
        title = "GO enrichment across clusters"
      ) +
      theme_bw(base_size = 9) +
      theme(
        axis.text.x = element_text(angle = 55, hjust = 1),
        panel.grid = element_blank()
      )
  } else {
    panel_b <- ggplot() +
      annotate(
        "text", x = 0, y = 0,
        label = paste0("No GO terms passed p < ", format(go_pvalue_cutoff, scientific = TRUE))
      ) +
      xlim(-1, 1) + ylim(-1, 1) + theme_void()
  }
} else {
  warning(
    "GO map not found: ", go_map_file,
    ". AgriGO-ready gene lists were exported, but enrichment was not calculated."
  )
  panel_b <- ggplot() +
    annotate(
      "text", x = 0, y = 0,
      label = "GO map unavailable. Use the exported AgriGO gene lists."
    ) +
    xlim(-1, 1) + ylim(-1, 1) + theme_void()
}

ggsave(
  file.path(figure_dir, "Figure_9B_GO_enrichment_heatmap.png"),
  panel_b, width = 13.5, height = 5.8, dpi = 300
)
ggsave(
  file.path(figure_dir, "Figure_9B_GO_enrichment_heatmap.pdf"),
  panel_b, width = 13.5, height = 5.8
)

figure_9 <- panel_a / panel_b +
  plot_annotation(tag_levels = "A") +
  plot_layout(heights = c(1.15, 1))

ggsave(
  file.path(figure_dir, "Figure_9_developmental_trends_and_GO.png"),
  figure_9, width = 13.5, height = 11.5, dpi = 300
)
ggsave(
  file.path(figure_dir, "Figure_9_developmental_trends_and_GO.pdf"),
  figure_9, width = 13.5, height = 11.5
)

# -----------------------------
# 7. Reproducibility record
# -----------------------------

writeLines(
  capture.output(sessionInfo()),
  file.path(session_dir, "07_developmental_trends_GO_sessionInfo.txt")
)

message("Developmental trend and GO workflow completed.")
message("Tables:  ", table_dir)
message("Figures: ", figure_dir)
