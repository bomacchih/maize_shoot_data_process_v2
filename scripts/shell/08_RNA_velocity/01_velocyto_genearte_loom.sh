#!/usr/bin/bash

# Define base paths using environment/user variables or relative paths

#XGE20_UL01_B73V5S210_short_MTCL/
#├── outs/
#│   ├── possorted_genome_bam.bam   # MUST be sorted & indexed
#│   ├── filtered_feature_bc_matrix/ (or raw_feature_bc_matrix/)
#│   └── barcodes.tsv.gz


#Expected Output
#Upon completion, velocyto will generate a .loom file inside the velocyto/ subfolder within each sample folder:
#XGE20122_short_MTCL/XGE20_UL01_B73V5S210_short_MTCL/velocyto/XGE20_UL01_B73V5S210_short_MTCL.loom

#!/usr/bin/env bash

set -euo pipefail  # Exit on error, unhandled variables, or pipe failures

# Define base directory paths
BASE_DIR="${HOME}/ST/spaceranger/spaceranger-2.1.0"
DATA_DIR="${BASE_DIR}/XGE202122_short_MTCL"
GTF_FILE="${BASE_DIR}/Zea_mays.Zm-B73-REFERENCE_NAM-5.0.51.gtf"

# Optional: Add repeat masker file path if available
# REPEAT_MASK="${BASE_DIR}/Zea_mays_rm.gtf"

# Define the sample IDs
SAMPLES=(
    "XGE20_UL01_B73V5S210_short_MTCL"
    "XGE20_UL02_B73V5S210_short_MTCL"
    "XGE20_UL03_B73V5S210_short_MTCL"
    "XGE20_UL04_B73V5S210_short_MTCL"
    "XGE21_VR01_B73V5S210_short_MTCL"
    "XGE21_VR02_B73V5S210_short_MTCL"
    "XGE21_VR03_B73V5S210_short_MTCL"
    "XGE21_VR04_B73V5S210_short_MTCL"
    "XGE22_DQ01_B73V5S210_short_MTCL"
    "XGE22_DQ02_B73V5S210_short_MTCL"
    "XGE22_DQ03_B73V5S210_short_MTCL"
    "XGE22_DQ04_B73V5S210_short_MTCL"
    "XGE22_DQ05_B73V5S210_short_MTCL"
    "XGE22_DQ06_B73V5S210_short_MTCL"
    "XGE22_DQ07_B73V5S210_short_MTCL"
    "XGE22_DQ08_B73V5S210_short_MTCL"
)

# Check if GTF file exists before starting
if [[ ! -f "${GTF_FILE}" ]]; then
    echo "Error: Reference GTF file not found at ${GTF_FILE}" >&2
    exit 1
fi

# Iterate over each sample
for SAMPLE in "${SAMPLES[@]}"; do
    SAMPLE_PATH="${DATA_DIR}/${SAMPLE}"
    BAM_FILE="${SAMPLE_PATH}/outs/possorted_genome_bam.bam"

    echo "=== Processing sample: ${SAMPLE} ==="

    # Check if sample directory exists
    if [[ ! -d "${SAMPLE_PATH}" ]]; then
        echo "Warning: Directory ${SAMPLE_PATH} does not exist. Skipping..."
        continue
    fi

    # Ensure BAM index exists before running velocyto
    if [[ -f "${BAM_FILE}" && ! -f "${BAM_FILE}.bai" ]]; then
        echo "Index missing for ${SAMPLE}. Indexing BAM file..."
        samtools index "${BAM_FILE}"
    fi

    # Run velocyto
    # Note: Append '-m ${REPEAT_MASK}' below if using a repeat masker GTF
    velocyto run10x "${SAMPLE_PATH}" "${GTF_FILE}"

    echo "Finished sample: ${SAMPLE}"
    echo "------------------------------------------------"
done