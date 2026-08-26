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
    # Choose a Seurat object currently loaded in the RStudio Global Environment.
    loaded_objects <- ls(envir = .GlobalEnv)
    seurat_objects <- loaded_objects[
        vapply(
            loaded_objects,
            function(object_name) {
                inherits(get(object_name, envir = .GlobalEnv), "Seurat")
            },
            logical(1)
        )
    ]

    seurat_object_name <- select.list(
        seurat_objects,
        multiple = FALSE,
        title = "Choose the loaded Seurat dataset"
    )

    # Choose the metadata CSV file interactively in RStudio.
    metadata_file <- file.choose()

    spot_metadata <- read.csv(
        metadata_file,
        stringsAsFactors = FALSE,
        check.names = FALSE
    )

    # Read the metadata type from the second column name.
    metadata_type <- colnames(spot_metadata)[2]

    # Match CSV barcodes to the spot order in the selected Seurat object.
    cell_ids <- colnames(get(seurat_object_name, envir = .GlobalEnv))
    metadata_values <- spot_metadata[[2]][
        match(cell_ids, spot_metadata$Barcode)
    ]
    names(metadata_values) <- cell_ids

    # Update the selected Seurat object directly under its existing object name.
    assignment <- substitute(
        TARGET[[METADATA_TYPE]] <- METADATA_VALUES,
        list(
            TARGET = as.name(seurat_object_name),
            METADATA_TYPE = metadata_type,
            METADATA_VALUES = metadata_values
        )
    )
    eval(assignment, envir = .GlobalEnv)

    # Show the imported metadata values.
    print(table(metadata_values, useNA = "ifany"))
})

