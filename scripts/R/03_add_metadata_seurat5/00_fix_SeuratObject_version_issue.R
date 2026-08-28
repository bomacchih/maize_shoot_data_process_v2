#!/usr/bin/env Rscript
# The Seurat RDS file downloaded from Zenodo uses an older SeuratObject
# data structure in which the spatial image lacks the "misc" slot.
# Run this update before importing metadata; otherwise, metadata transfer
# and removal of spots with missing metadata may fail.


library(Seurat)
library(SeuratObject)

# Keep temporary objects inside a local environment.
local({
    # Find Seurat objects loaded in the Global Environment.
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
        stop("No Seurat objects are loaded in the Global Environment.")
    }

    # Select the object to update.
    seurat_object_name <- select.list(
        seurat_objects,
        multiple = FALSE,
        title = "Choose the Seurat object to update"
    )

    if (!nzchar(seurat_object_name)) {
        stop("No Seurat object was selected.")
    }

    target_object <- get(
        seurat_object_name,
        envir = .GlobalEnv
    )

    # Preserve the active identities.
    idents_before <- Idents(target_object)

    # Update each spatial image object.
    image_names <- names(target_object@images)

    if (length(image_names) > 0) {
        for (image_name in image_names) {
            target_object@images[[image_name]] <-
                SeuratObject::UpdateSlots(
                    target_object@images[[image_name]]
                )
        }
    }

    # Update the complete Seurat object.
    target_object <-
        SeuratObject::UpdateSlots(target_object)

    target_object <-
        SeuratObject::UpdateSeuratObject(target_object)

    # Confirm that the cells have not changed.
    if (!setequal(names(idents_before), colnames(target_object))) {
        stop("Cell names changed unexpectedly during the update.")
    }

    # Confirm that active identity values have not changed.
    stopifnot(
        identical(
            as.character(
                Idents(target_object)[names(idents_before)]
            ),
            as.character(idents_before)
        )
    )

    # Validate the updated object.
    validObject(target_object)

    # Replace the selected Global Environment object.
    assign(
        seurat_object_name,
        target_object,
        envir = .GlobalEnv
    )

    message(
        "Successfully updated: ",
        seurat_object_name
    )

    message(
        "Stored object version: ",
        as.character(target_object@version)
    )

    if (length(image_names) > 0) {
        message(
            "Updated image objects: ",
            paste(image_names, collapse = ", ")
        )

        print(
            SpatialDimPlot(
                target_object,
                images = image_names[1]
            )
        )
    } else {
        message("The selected object does not contain spatial images.")
    }
})
