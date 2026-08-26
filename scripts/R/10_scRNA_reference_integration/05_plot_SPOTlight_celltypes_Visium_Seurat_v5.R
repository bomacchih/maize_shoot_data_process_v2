# Replot selected SPOTlight cell-type proportions without rerunning deconvolution
#
# Input:
#   an in-memory Seurat object named `visium_mapped`, or
#   data/processed/maize_shoot_14samples_celltype_mapped_SPOTlight_seurat_v5.rds

suppressPackageStartupMessages({
  library(Seurat)
  library(SPOTlight)
  library(ggplot2)
  library(patchwork)
})

vascular_cutoff <- 0.10
sam_cutoff <- 0.05

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

first_existing <- function(candidates, available, required = TRUE) {
  selected <- candidates[candidates %in% available]
  if (length(selected)) return(selected[[1L]])
  if (required) stop("None of these fields was found: ",
                     paste(candidates, collapse = ", "))
  NA_character_
}

project_root <- find_project_root()
input_rds <- file.path(
  project_root, "data", "processed",
  "maize_shoot_14samples_celltype_mapped_SPOTlight_seurat_v5.rds"
)
existing_subset_rds <- file.path(
  project_root, "data", "processed",
  "XGE202122_S5_subset_embleaf_harmony_join.rds"
)
figure_dir <- file.path(
  project_root, "results", "figures", "12_scRNA_Visium_mapping",
  "selected_celltypes"
)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

if (exists("visium_mapped", envir = .GlobalEnv, inherits = FALSE)) {
  object <- get("visium_mapped", envir = .GlobalEnv)
} else if (file.exists(input_rds)) {
  object <- readRDS(input_rds)
} else if (file.exists(existing_subset_rds)) {
  object <- readRDS(existing_subset_rds)
} else {
  stop("No mapped Visium object was found.")
}
stopifnot(inherits(object, "Seurat"))
original_idents <- Idents(object)
metadata_columns <- colnames(object[[]])

vascular_feature <- first_existing(
  c("SPOT_Vascular_tissue", "SPOT_Vascular"),
  metadata_columns, required = FALSE
)
sam_feature <- first_existing(
  c("SPOT_SAM_tissue", "SPOT_Shoot_apical_meristem", "SPOT_SAM"),
  metadata_columns, required = FALSE
)
selected_features <- unique(na.omit(c(
  vascular_feature,
  sam_feature,
  intersect(c("SPOT_Mesophyll", "SPOT_Bundle_sheath"), metadata_columns)
)))
if (!length(selected_features)) {
  stop("No vascular or SAM SPOTlight proportion fields were found.")
}

umap_reduction <- first_existing(
  c("umapharmony", "umap.harmony", "harmony.umap", "umap", "ref.umap"),
  Reductions(object)
)

umap_proportions <- FeaturePlot(
  object,
  features = selected_features,
  reduction = umap_reduction,
  order = TRUE,
  min.cutoff = 0,
  max.cutoff = 1,
  ncol = length(selected_features)
) + plot_annotation(title = "Selected SPOTlight proportions on Visium UMAP")

ggsave(
  file.path(figure_dir, "selected_SPOTlight_proportions_UMAP.png"),
  umap_proportions, width = 6 * length(selected_features), height = 5.5,
  dpi = 300
)
ggsave(
  file.path(figure_dir, "selected_SPOTlight_proportions_UMAP.pdf"),
  umap_proportions, width = 6 * length(selected_features), height = 5.5
)

threshold_plot <- function(object, feature, cutoff, title) {
  high_cells <- colnames(object)[
    is.finite(object[[feature, drop = TRUE]]) &
      object[[feature, drop = TRUE]] > cutoff
  ]
  DimPlot(
    object,
    reduction = umap_reduction,
    cells.highlight = high_cells,
    cols = "grey82",
    cols.highlight = "#E41A1C",
    pt.size = 0.25,
    sizes.highlight = 0.45
  ) +
    NoLegend() +
    ggtitle(paste0(title, "\n", feature, " > ", cutoff)) +
    theme(plot.title = element_text(size = 12, hjust = 0.5, face = "bold"))
}

threshold_plots <- list()
if (!is.na(vascular_feature)) {
  threshold_plots[["vascular"]] <- threshold_plot(
    object, vascular_feature, vascular_cutoff, "Vascular-enriched spots"
  )
}
if (!is.na(sam_feature)) {
  threshold_plots[["SAM"]] <- threshold_plot(
    object, sam_feature, sam_cutoff, "SAM-enriched spots"
  )
}
threshold_figure <- wrap_plots(threshold_plots, ncol = 2)
ggsave(
  file.path(figure_dir, "selected_SPOTlight_thresholds_UMAP.png"),
  threshold_figure, width = 11, height = 5.5, dpi = 300
)
ggsave(
  file.path(figure_dir, "selected_SPOTlight_thresholds_UMAP.pdf"),
  threshold_figure, width = 11, height = 5.5
)

# Figure C: continuous SPOTlight proportions for all 12 cell types.
figure_c_celltypes <- c(
  "Vascular_tissue", "Leaf_rim", "Shoot_system_epidermis",
  "Leaf_primordium", "Pavement_cell_N", "Leaf_guard_cell",
  "Bundle_sheath", "Mesophyll", "Leaf_epidermis",
  "Leaf_subsidiary_cell", "Shoot_apical_meristem", "Pavement_cell_A"
)
figure_c_features <- paste0("SPOT_", figure_c_celltypes)
figure_c_features <- figure_c_features[figure_c_features %in% metadata_columns]
umap_coordinates <- as.data.frame(Embeddings(object, reduction = umap_reduction))
colnames(umap_coordinates)[1:2] <- c("UMAP_1", "UMAP_2")

figure_c_plots <- lapply(figure_c_features, function(feature) {
  plot_data <- data.frame(
    UMAP_1 = umap_coordinates$UMAP_1,
    UMAP_2 = umap_coordinates$UMAP_2,
    proportion = object[[feature, drop = TRUE]]
  )
  plot_data$proportion[!is.finite(plot_data$proportion)] <- 0
  plot_data <- plot_data[order(plot_data$proportion), , drop = FALSE]
  ggplot(plot_data, aes(UMAP_1, UMAP_2, color = proportion)) +
    geom_point(size = 0.28, alpha = 0.9) +
    scale_color_gradientn(
      colors = c("#E5E5E5", "#FCAE91", "#FB6A4A", "#CB181D"),
      limits = c(0, 1),
      oob = scales::squish
    ) +
    coord_equal() +
    labs(
      title = sub("^SPOT_", "", feature),
      x = "umapharmony_1",
      y = "umapharmony_2"
    ) +
    theme_classic(base_size = 8) +
    theme(
      legend.position = "none",
      plot.title = element_text(face = "bold", hjust = 0.5, size = 8),
      axis.title = element_text(size = 8),
      axis.text = element_text(size = 7)
    )
})

figure_c <- wrap_plots(figure_c_plots, ncol = 4) +
  plot_annotation(
    title = "C   SPOTlight-estimated cell-type proportions on Visium UMAP"
  )
ggsave(
  file.path(figure_dir, "Figure_C_SPOTlight_proportion_UMAPs.png"),
  figure_c, width = 12, height = 9, dpi = 300
)
ggsave(
  file.path(figure_dir, "Figure_C_SPOTlight_proportion_UMAPs.pdf"),
  figure_c, width = 12, height = 9
)

# Spatial feature maps are produced per image so spots are overlaid on the
# appropriate capture-area image. Each PDF contains one page per image.
image_names <- Images(object)
if (length(image_names)) {
  for (feature in selected_features) {
    pdf(
      file.path(figure_dir, paste0(feature, "_spatial_maps_by_image.pdf")),
      width = 7.5, height = 7.5, useDingbats = FALSE
    )
    for (image_name in image_names) {
      print(
        SpatialFeaturePlot(
          object,
          features = feature,
          images = image_name,
          min.cutoff = 0,
          max.cutoff = 1,
          pt.size.factor = 1.6
        ) + ggtitle(paste(feature, image_name, sep = ": "))
      )
    }
    dev.off()
  }

  # Representative spatial scatter-pie from an existing capture area.
  representative_image <- if ("VR03" %in% image_names) "VR03" else image_names[[1L]]
  coordinates <- GetTissueCoordinates(object, image = representative_image)
  representative_spots <- intersect(rownames(coordinates), colnames(object))
  pie_features <- grep("^SPOT_", metadata_columns, value = TRUE)
  pie_features <- setdiff(
    pie_features,
    c("SPOT_top_type", "SPOT_top_prop", "SPOT_high_purity",
      "SPOT_entropy", "SPOT_entropy_normalized")
  )
  if (length(representative_spots) && length(pie_features)) {
    proportions <- as.matrix(
      object[[]][representative_spots, pie_features, drop = FALSE]
    )
    colnames(proportions) <- sub("^SPOT_", "", colnames(proportions))
    proportions[!is.finite(proportions) | proportions < 0] <- 0
    row_totals <- rowSums(proportions)
    valid_spots <- row_totals > 0
    proportions <- proportions[valid_spots, , drop = FALSE]
    proportions <- proportions / rowSums(proportions)
    coordinates <- coordinates[rownames(proportions), , drop = FALSE]

    coordinate_pairs <- list(
      c("imagecol", "imagerow"), c("x", "y"), c("col", "row"),
      c("pxl_col_in_fullres", "pxl_row_in_fullres")
    )
    xy_columns <- NULL
    for (pair in coordinate_pairs) {
      if (all(pair %in% colnames(coordinates))) {
        xy_columns <- pair
        break
      }
    }
    if (is.null(xy_columns)) {
      numeric_columns <- colnames(coordinates)[
        vapply(coordinates, is.numeric, logical(1))
      ]
      xy_columns <- numeric_columns[1:2]
    }

    scatterpie_plot <- SPOTlight::plotSpatialScatterpie(
      x = as.matrix(coordinates[, xy_columns, drop = FALSE]),
      y = proportions,
      cell_types = colnames(proportions),
      img = FALSE,
      scatterpie_alpha = 0.85,
      pie_scale = 0.45
    ) +
      ggtitle(paste("SPOTlight proportions:", representative_image)) +
      theme(
        plot.background = element_rect(fill = "white", color = NA),
        panel.background = element_rect(fill = "white", color = NA),
        plot.title = element_text(color = "black", face = "bold", hjust = 0.5),
        legend.text = element_text(color = "black"),
        legend.title = element_text(color = "black")
      )
    ggsave(
      file.path(figure_dir, "representative_VR03_SPOTlight_scatterpie.png"),
      scatterpie_plot, width = 9, height = 7, dpi = 300
    )
    ggsave(
      file.path(figure_dir, "representative_VR03_SPOTlight_scatterpie.pdf"),
      scatterpie_plot, width = 9, height = 7
    )

    # Figure B: the four VR03 sections plotted separately with a shared legend.
    section_column <- first_existing(
      c("section", "section_id"), metadata_columns, required = FALSE
    )
    if (!is.na(section_column)) {
      section_values <- as.character(
        object[[]][representative_spots, section_column, drop = TRUE]
      )
      names(section_values) <- representative_spots
      section_levels <- unique(section_values[!is.na(section_values) &
                                                nzchar(section_values)])
      section_number <- suppressWarnings(as.numeric(
        sub(".*?([0-9]+)$", "\\1", section_levels)
      ))
      section_levels <- section_levels[order(section_number, section_levels)]
      section_levels <- head(section_levels, 4L)

      figure_b_celltypes <- c(
        "Bundle_sheath", "Leaf_epidermis", "Leaf_guard_cell",
        "Leaf_primordium", "Leaf_rim", "Mesophyll", "Pavement_cell_A",
        "Shoot_apical_meristem", "Shoot_system_epidermis"
      )
      figure_b_features <- paste0("SPOT_", figure_b_celltypes)
      figure_b_features <- figure_b_features[
        figure_b_features %in% metadata_columns
      ]
      figure_b_celltypes <- sub("^SPOT_", "", figure_b_features)
      figure_b_palette <- c(
        Bundle_sheath = "#1B1B1B",
        Leaf_epidermis = "#1B7F79",
        Leaf_guard_cell = "#8C6BB1",
        Leaf_primordium = "#E78AC3",
        Leaf_rim = "#5E3C99",
        Mesophyll = "#9ECAE1",
        Pavement_cell_A = "#377EB8",
        Shoot_apical_meristem = "#B15928",
        Shoot_system_epidermis = "#33A02C"
      )

      section_plots <- lapply(seq_along(section_levels), function(index) {
        section_name <- section_levels[[index]]
        section_spots <- names(section_values)[section_values == section_name]
        section_spots <- intersect(section_spots, rownames(coordinates))
        section_proportions <- as.matrix(
          object[[]][section_spots, figure_b_features, drop = FALSE]
        )
        colnames(section_proportions) <- figure_b_celltypes
        section_proportions[
          !is.finite(section_proportions) | section_proportions < 0
        ] <- 0
        keep <- rowSums(section_proportions) > 0
        section_proportions <- section_proportions[keep, , drop = FALSE]
        section_proportions <- section_proportions /
          rowSums(section_proportions)
        section_coordinates <- coordinates[
          rownames(section_proportions), xy_columns, drop = FALSE
        ]

        SPOTlight::plotSpatialScatterpie(
          x = as.matrix(section_coordinates),
          y = section_proportions,
          cell_types = figure_b_celltypes,
          img = FALSE,
          scatterpie_alpha = 0.9,
          pie_scale = 0.75
        ) +
          scale_fill_manual(
            values = figure_b_palette,
            breaks = figure_b_celltypes,
            drop = FALSE
          ) +
          ggtitle(paste0("S", index)) +
          theme(
            plot.background = element_rect(fill = "white", color = NA),
            panel.background = element_rect(fill = "white", color = NA),
            plot.title = element_text(
              color = "black", face = "bold", hjust = 0.5, size = 13
            ),
            legend.position = "none"
          )
      })

      # Build one explicit legend because plotSpatialScatterpie can emit
      # duplicate fill guides when several scatter-pie panels are combined.
      legend_data <- data.frame(
        cell_type = factor(figure_b_celltypes, levels = figure_b_celltypes),
        y = rev(seq_along(figure_b_celltypes))
      )
      legend_panel <- ggplot(legend_data, aes(x = 0, y = y)) +
        geom_tile(
          aes(fill = cell_type), width = 0.18, height = 0.58,
          color = NA
        ) +
        geom_text(
          aes(x = 0.16, label = cell_type), hjust = 0,
          size = 3.2, color = "black"
        ) +
        scale_fill_manual(values = figure_b_palette, guide = "none") +
        coord_cartesian(xlim = c(-0.12, 1.85), clip = "off") +
        labs(title = "Cell type") +
        theme_void() +
        theme(
          plot.title = element_text(face = "bold", hjust = 0, size = 10),
          plot.margin = margin(20, 5, 5, 5)
        )

      figure_b <- wrap_plots(
        c(section_plots, list(legend_panel)),
        nrow = 1,
        widths = c(rep(1, length(section_plots)), 0.95)
      ) +
        plot_annotation(
          title = "B   SPOTlight deconvolution across VR03 sections"
        )
      ggsave(
        file.path(figure_dir, "Figure_B_VR03_section_scatterpies.png"),
        figure_b, width = 14, height = 4.6, dpi = 300
      )
      ggsave(
        file.path(figure_dir, "Figure_B_VR03_section_scatterpies.pdf"),
        figure_b, width = 14, height = 4.6
      )
    }
  }
}

Idents(object) <- original_idents
stopifnot(identical(Idents(object), original_idents))
message("Selected SPOTlight plots completed: ", figure_dir)
