#!/usr/bin/env Rscript

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
# │       └── UL01_stutility.rds
# └── scripts/
#     └── R/
#         └── 02_load_spaceranger/
#             └── 01_load_single_sample_stutility.R
#
# Replace UL01 with the sample ID being processed, such as VR01 or DQ08.

library(Seurat)
library(STutility)

# Manually enter one sample ID, for example UL01, VR01, or DQ08.
sample_id <- "UL01"

# Space Ranger output directory for this capture area.
sample_dir <- file.path("data", "raw", sample_id, "outs")

samples <- file.path(sample_dir, "filtered_feature_bc_matrix.h5")
imgs <- file.path(sample_dir, "spatial", "tissue_hires_image.png")
spotfiles <- file.path(sample_dir, "spatial", "tissue_positions.csv")
json <- file.path(sample_dir, "spatial", "scalefactors_json.json")

infoTable <- data.frame(
    samples = samples,
    imgs = imgs,
    spotfiles = spotfiles,
    json = json,
    sample_id = sample_id,
    stringsAsFactors = FALSE
)

se <- InputFromTable(
    infotable = infoTable,
    platform = "Visium"
)

dir.create(file.path("data", "processed"), recursive = TRUE, showWarnings = FALSE)

saveRDS(
    se,
    file = file.path("data", "processed", paste0(sample_id, "_stutility.rds"))
)


session_dir <- file.path("results", "sessionInfo")
dir.create(session_dir, recursive = TRUE, showWarnings = FALSE)
writeLines(
    capture.output(sessionInfo()),
    file.path(session_dir, "02_load_spaceranger_sessionInfo.txt")
)

