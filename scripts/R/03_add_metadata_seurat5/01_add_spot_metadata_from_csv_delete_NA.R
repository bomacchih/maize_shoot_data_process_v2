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
# Keep temporary objects inside a local environment.
# Require new imported seurat data from raw files. rds downloaded from Zenodo has image misc info issue.

local({

    # Find loaded Seurat objects.
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
        stop("No Seurat object was found in the Global Environment.")
    }

    # Select the Seurat object.
    seurat_object_name <- select.list(
        seurat_objects,
        multiple = FALSE,
        title = "Choose the loaded Seurat dataset"
    )

    if (!nzchar(seurat_object_name)) {
        stop("No Seurat object was selected.")
    }

    seurat_object <- get(
        seurat_object_name,
        envir = .GlobalEnv
    )

    # Select and read the metadata CSV file.
    metadata_file <- file.choose()

    spot_metadata <- read.csv(
        metadata_file,
        stringsAsFactors = FALSE,
        check.names = FALSE
    )

    if (!"Barcode" %in% colnames(spot_metadata)) {
        stop("The CSV file does not contain a Barcode column.")
    }

    # Metadata columns to import.
    columns_to_import <- c(
        "orig.ident",
        "section_id",
        "domains",
        "sample_id",
        "section",
        "cca_clusters",
        "seurat_clusters",
        "sample_domain",
        "ms_ve",
        "domain_section"
    )

    missing_columns <- setdiff(
        columns_to_import,
        colnames(spot_metadata)
    )

    if (length(missing_columns) > 0) {
        stop(
            paste(
                "Columns missing from the CSV:",
                paste(missing_columns, collapse = ", ")
            )
        )
    }

    # Remove suffixes such as _1_6 from barcodes.
    normalize_barcode <- function(x) {
        sub("_[0-9]+_[0-9]+$", "", as.character(x))
    }

    cell_ids <- colnames(seurat_object)

    # Test several possible barcode-matching methods.
    matching_results <- list(
        Barcode_exact = match(
            cell_ids,
            spot_metadata$Barcode
        ),

        Barcode_normalized = match(
            normalize_barcode(cell_ids),
            normalize_barcode(spot_metadata$Barcode)
        )
    )

    # spot_inS5 may retain the original barcode suffix.
    if ("spot_inS5" %in% colnames(spot_metadata)) {
        matching_results$spot_inS5_exact <- match(
            cell_ids,
            spot_metadata$spot_inS5
        )

        matching_results$spot_inS5_normalized <- match(
            normalize_barcode(cell_ids),
            normalize_barcode(spot_metadata$spot_inS5)
        )
    }

    # Count matches under each method.
    match_counts <- vapply(
        matching_results,
        function(x) sum(!is.na(x)),
        integer(1)
    )

    print(match_counts)

    # Use the method giving the largest number of matches.
    selected_method <- names(which.max(match_counts))
    metadata_index <- matching_results[[selected_method]]

    message("Matching method used: ", selected_method)
    message("Matched spots: ", sum(!is.na(metadata_index)))
    message("Unmatched spots: ", sum(is.na(metadata_index)))

    if (sum(!is.na(metadata_index)) == 0) {
        stop("No Seurat spots matched the metadata file.")
    }

    # Identify unmatched Seurat spots.
    matched_spots <- !is.na(metadata_index)

    unmatched_cell_ids <- cell_ids[!matched_spots]
    matched_cell_ids <- cell_ids[matched_spots]

    if (length(unmatched_cell_ids) > 0) {
        message(
            "Removing ",
            length(unmatched_cell_ids),
            " unmatched spots from the Seurat object."
        )

        message("First unmatched spot IDs:")

        print(
            head(unmatched_cell_ids, 20)
        )
    }

    # Keep only spots present in the metadata file.
    seurat_object <- subset(
        seurat_object,
        cells = matched_cell_ids
    )

    metadata_index <- metadata_index[matched_spots]
    cell_ids <- matched_cell_ids

    # Confirm that spot order remains correct.
    if (!identical(colnames(seurat_object), cell_ids)) {
        stop("Spot order changed unexpectedly after subsetting.")
    }

    # Create the metadata table in Seurat spot order.
    metadata_to_add <- spot_metadata[
        metadata_index,
        columns_to_import,
        drop = FALSE
    ]

    rownames(metadata_to_add) <- cell_ids

    # Store cluster labels as categorical variables.
    cluster_columns <- intersect(
        c("cca_clusters", "seurat_clusters"),
        colnames(metadata_to_add)
    )

    metadata_to_add[cluster_columns] <- lapply(
        metadata_to_add[cluster_columns],
        factor
    )

    # Add metadata to the Seurat object.
    seurat_object <- SeuratObject::AddMetaData(
        object = seurat_object,
        metadata = metadata_to_add
    )

    # Replace the original object in the Global Environment.
    assign(
        seurat_object_name,
        seurat_object,
        envir = .GlobalEnv
    )

    # Final report.
    message("Updated Seurat object: ", seurat_object_name)
    message("Remaining spots: ", ncol(seurat_object))
    message("Removed unmatched spots: ", length(unmatched_cell_ids))

    cat("\nNA count for each imported column:\n")

    print(
        colSums(
            is.na(seurat_object@meta.data[, columns_to_import, drop = FALSE])
        )
    )
})


session_dir <- file.path("results", "sessionInfo")
dir.create(session_dir, recursive = TRUE, showWarnings = FALSE)
writeLines(
    capture.output(sessionInfo()),
    file.path(session_dir, "03_add_metadata_seurat5_sessionInfo.txt")
)
