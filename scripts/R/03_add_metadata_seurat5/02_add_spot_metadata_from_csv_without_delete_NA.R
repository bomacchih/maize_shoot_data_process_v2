#!/usr/bin/env Rscript

# Add spot metadata from a CSV file directly to a Seurat v5 object in RStudio.
# The selected Seurat object is updated under its existing name. No duplicate
# Seurat object or output file is created. Temporary objects are removed when
# the script finishes.
#
# Before running this script, manually load one or more Seurat objects:
#
# UL01 <- readRDS("data/processed/UL01_seurat_v5.rds")
# VR01 <- readRDS("data/processed/VR01_seurat_v5.rds")
#
# Expected metadata file:
# data/metadata/UL01_metadata/filename.csv
#
# The first CSV column must be named Barcode. The second column contains the
# metadata values, and its column name will be used as the metadata type.

library(Seurat)

# Keep all temporary objects inside a local environment.
local({
    # Identify Seurat objects in the Global Environment.
    loaded_objects <- ls(envir = .GlobalEnv)

    seurat_objects <- loaded_objects[
        vapply(
            loaded_objects,
            function(object_name) {
                inherits(
                    get(object_name, envir = .GlobalEnv),
                    "Seurat"
                )
            },
            logical(1)
        )
    ]

    if (length(seurat_objects) == 0) {
        stop("No Seurat object is loaded in the Global Environment.")
    }

    # Choose the Seurat object.
    seurat_object_name <- select.list(
        seurat_objects,
        multiple = FALSE,
        title = "Choose the loaded Seurat dataset"
    )

    if (!nzchar(seurat_object_name)) {
        stop("No Seurat object was selected.")
    }

    target_object <- get(
        seurat_object_name,
        envir = .GlobalEnv
    )

    # Preserve the active identities.
    active_ident_before <- SeuratObject::Idents(target_object)

    # Choose and read the metadata CSV file.
    metadata_file <- file.choose()

    spot_metadata <- read.csv(
        metadata_file,
        stringsAsFactors = FALSE,
        check.names = FALSE
    )

    if (!"Barcode" %in% colnames(spot_metadata)) {
        stop("The metadata file must contain a column named 'Barcode'.")
    }

    spot_metadata$Barcode <- trimws(
        as.character(spot_metadata$Barcode)
    )

    if (anyDuplicated(spot_metadata$Barcode)) {
        stop("Duplicated barcodes were found in the metadata CSV.")
    }

    # Extract the original 10x barcode from Seurat cell names.
    # Handles:
    # AAACAAGTATCTCCCA-1
    # VR03_AAACAAGTATCTCCCA-1
    # AAACAAGTATCTCCCA-1_1_6
    extract_10x_barcode <- function(cell_names) {
        barcode_pattern <- "[ACGTN]+-[0-9]+"

        has_barcode <- grepl(
            barcode_pattern,
            cell_names,
            perl = TRUE
        )

        extracted <- rep(NA_character_, length(cell_names))

        extracted[has_barcode] <- sub(
            paste0("^.*?(", barcode_pattern, ").*$"),
            "\\1",
            cell_names[has_barcode],
            perl = TRUE
        )

        extracted
    }

    cell_ids <- colnames(target_object)
    cell_barcodes <- extract_10x_barcode(cell_ids)

    # Restrict matching to the sample represented by the CSV when
    # the selected object contains multiple samples.
    target_cells <- rep(TRUE, length(cell_ids))

    if (
        "sample_id" %in% colnames(spot_metadata) &&
        "sample_id" %in% colnames(target_object[[]])
    ) {
        csv_sample_ids <- unique(
            na.omit(as.character(spot_metadata$sample_id))
        )

        if (length(csv_sample_ids) == 1) {
            target_cells <- (
                as.character(target_object$sample_id) ==
                csv_sample_ids
            )
        }
    }

    # Prevent ambiguous matching across multiple samples.
    if (
        all(target_cells) &&
        anyDuplicated(cell_barcodes[!is.na(cell_barcodes)])
    ) {
        stop(
            paste(
                "Raw barcodes are duplicated in the selected Seurat object.",
                "Use an object containing one sample or ensure that",
                "meta.data$sample_id is available."
            )
        )
    }

    csv_index <- match(
        cell_barcodes,
        spot_metadata$Barcode
    )

    matched_cells <- target_cells & !is.na(csv_index)

    if (!any(matched_cells)) {
        stop(
            "No matching barcodes were found between the Seurat object and CSV."
        )
    }

    # Import every CSV column except Barcode.
    metadata_columns <- setdiff(
        colnames(spot_metadata),
        "Barcode"
    )

    for (metadata_type in metadata_columns) {
        source_values <- spot_metadata[[metadata_type]]

        if (is.factor(source_values)) {
            source_values <- as.character(source_values)
        }

        # Preserve existing values for unmatched cells.
        if (metadata_type %in% colnames(target_object[[]])) {
            metadata_values <- target_object[[]][[metadata_type]]

            if (is.factor(metadata_values)) {
                metadata_values <- as.character(metadata_values)
            }
        } else {
            metadata_values <- rep(
                NA,
                length(cell_ids)
            )
        }

        metadata_values[matched_cells] <-
            source_values[csv_index[matched_cells]]

        names(metadata_values) <- cell_ids

        target_object[[metadata_type]] <- metadata_values
    }

    # Confirm that metadata import did not change active identities.
    if (
        !identical(
            as.character(SeuratObject::Idents(target_object)),
            as.character(active_ident_before)
        )
    ) {
        stop("The active identities changed unexpectedly.")
    }

    # Replace the selected object without creating another global object.
    assign(
        seurat_object_name,
        target_object,
        envir = .GlobalEnv
    )

    message(
        sum(matched_cells),
        " of ",
        sum(target_cells),
        " target spots matched the CSV."
    )

    message(
        "Imported metadata columns: ",
        paste(metadata_columns, collapse = ", ")
    )
})


session_dir <- file.path("results", "sessionInfo")
dir.create(session_dir, recursive = TRUE, showWarnings = FALSE)
writeLines(
    capture.output(sessionInfo()),
    file.path(session_dir, "03_add_metadata_seurat5_sessionInfo.txt")
)
