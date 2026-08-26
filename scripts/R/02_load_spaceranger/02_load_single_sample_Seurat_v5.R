#!/usr/bin/env Rscript

# Load one Space Ranger capture-area output with Seurat v5.
#
# Expected file structure (run this script from the repository root):
#
# maize_shoot_data_process_v2/
# ├── data/
# │   ├── raw/
# │   │   └── UL01/
# │   │       └── outs/
# │   │           ├── filtered_feature_bc_matrix.h5
# │   │           └── spatial/
# │   │               ├── tissue_hires_image.png
# │   │               ├── tissue_positions_list.csv
# │   │               └── scalefactors_json.json
# │   └── processed/
# │       └── UL01_seurat_v5.rds
# └── scripts/
#     └── R/
#         └── 02_load_spaceranger/
#             └── 02_load_single_sample_Seurat_v5.R

library(Seurat)

# Manually enter one sample ID, for example UL01, VR01, or DQ08.
sample_id <- "UL01"

# Space Ranger output directory for this capture area.
data_dir <- file.path("data", "raw", sample_id, "outs")

se <- Load10X_Spatial(
    data.dir = data_dir,
    filename = "filtered_feature_bc_matrix.h5",
    assay = "RNA",
    slice = sample_id,
    image.name = "tissue_hires_image.png"
)

se$sample_id <- sample_id
se$orig.ident <- sample_id

dir.create(file.path("data", "processed"), recursive = TRUE, showWarnings = FALSE)

saveRDS(
    se,
    file = file.path("data", "processed", paste0(sample_id, "_seurat_v5.rds"))
)

