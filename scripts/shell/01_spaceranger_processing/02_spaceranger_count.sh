#!/usr/bin/env bash
set -euo pipefail

# Run Space Ranger count for one Visium capture area.
#
# This script contains no experiment-specific absolute paths. Relative input
# paths are resolved from the repository root, so the script can be called from
# any working directory. Run the script once for each capture area.
#
# Example:
# bash scripts/shell/01_spaceranger_processing/02_spaceranger_count.sh \
#   --id XGE21_VR03_B73V5 \
#   --sample XGE21-VR03 \
#   --transcriptome data/reference/B73_V5 \
#   --fastqs data/raw/fastq/XGE21-VR03 \
#   --image data/raw/images/XGE21-VR03.tif \
#   --slide V10F04-109 \
#   --area C1 \
#   --output-dir data/raw/spaceranger \
#   --localcores 8 \
#   --localmem 64

usage() {
    cat <<'EOF'
Usage:
  02_spaceranger_count.sh \
    --id RUN_ID \
    --sample FASTQ_SAMPLE_NAME \
    --transcriptome REFERENCE_DIRECTORY \
    --fastqs FASTQ_DIRECTORY \
    --image HISTOLOGY_IMAGE \
    --slide SLIDE_ID \
    --area CAPTURE_AREA \
    [--output-dir OUTPUT_DIRECTORY] \
    [--localcores INTEGER] \
    [--localmem INTEGER]

Required arguments:
  --id              Name of the Space Ranger output directory.
  --sample          Sample prefix used in the FASTQ filenames.
  --transcriptome   Space Ranger reference directory.
  --fastqs          Directory containing the FASTQ files.
  --image           Brightfield histology image file.
  --slide           Visium slide serial number.
  --area            Visium capture area, for example A1 or C1.

Optional arguments:
  --output-dir      Parent output directory.
                    Default: data/raw/spaceranger
  --localcores      Number of CPU cores. Default: 8
  --localmem        Memory in GB. Default: 64
  --help            Show this message.

Environment variable:
  SPACERANGER_BIN   Space Ranger executable or command name.
                    Default: spaceranger
EOF
}

require_argument_value() {
    local option_name="$1"
    local option_value="${2:-}"

    if [[ -z "${option_value}" || "${option_value}" == --* ]]; then
        echo "Error: ${option_name} requires a value." >&2
        usage >&2
        exit 2
    fi
}

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd -- "${SCRIPT_DIR}/../../.." && pwd)"

SPACERANGER_BIN="${SPACERANGER_BIN:-spaceranger}"
OUTPUT_DIR="data/raw/spaceranger"
LOCALCORES=8
LOCALMEM=64

RUN_ID=""
SAMPLE_NAME=""
TRANSCRIPTOME=""
FASTQ_DIR=""
IMAGE_FILE=""
SLIDE_ID=""
CAPTURE_AREA=""

while [[ $# -gt 0 ]]; do
    case "$1" in
        --id)
            require_argument_value "$1" "${2:-}"
            RUN_ID="$2"
            shift 2
            ;;
        --sample)
            require_argument_value "$1" "${2:-}"
            SAMPLE_NAME="$2"
            shift 2
            ;;
        --transcriptome)
            require_argument_value "$1" "${2:-}"
            TRANSCRIPTOME="$2"
            shift 2
            ;;
        --fastqs)
            require_argument_value "$1" "${2:-}"
            FASTQ_DIR="$2"
            shift 2
            ;;
        --image)
            require_argument_value "$1" "${2:-}"
            IMAGE_FILE="$2"
            shift 2
            ;;
        --slide)
            require_argument_value "$1" "${2:-}"
            SLIDE_ID="$2"
            shift 2
            ;;
        --area)
            require_argument_value "$1" "${2:-}"
            CAPTURE_AREA="$2"
            shift 2
            ;;
        --output-dir)
            require_argument_value "$1" "${2:-}"
            OUTPUT_DIR="$2"
            shift 2
            ;;
        --localcores)
            require_argument_value "$1" "${2:-}"
            LOCALCORES="$2"
            shift 2
            ;;
        --localmem)
            require_argument_value "$1" "${2:-}"
            LOCALMEM="$2"
            shift 2
            ;;
        --help|-h)
            usage
            exit 0
            ;;
        *)
            echo "Error: unknown argument: $1" >&2
            usage >&2
            exit 2
            ;;
    esac
done

missing_arguments=()
[[ -z "${RUN_ID}" ]] && missing_arguments+=("--id")
[[ -z "${SAMPLE_NAME}" ]] && missing_arguments+=("--sample")
[[ -z "${TRANSCRIPTOME}" ]] && missing_arguments+=("--transcriptome")
[[ -z "${FASTQ_DIR}" ]] && missing_arguments+=("--fastqs")
[[ -z "${IMAGE_FILE}" ]] && missing_arguments+=("--image")
[[ -z "${SLIDE_ID}" ]] && missing_arguments+=("--slide")
[[ -z "${CAPTURE_AREA}" ]] && missing_arguments+=("--area")

if [[ ${#missing_arguments[@]} -gt 0 ]]; then
    echo "Error: missing required arguments: ${missing_arguments[*]}" >&2
    usage >&2
    exit 2
fi

resolve_from_project_root() {
    local input_path="$1"

    if [[ "${input_path}" = /* ]]; then
        printf '%s\n' "${input_path}"
    else
        printf '%s/%s\n' "${PROJECT_ROOT}" "${input_path#./}"
    fi
}

TRANSCRIPTOME="$(resolve_from_project_root "${TRANSCRIPTOME}")"
FASTQ_DIR="$(resolve_from_project_root "${FASTQ_DIR}")"
IMAGE_FILE="$(resolve_from_project_root "${IMAGE_FILE}")"
OUTPUT_DIR="$(resolve_from_project_root "${OUTPUT_DIR}")"

if [[ ! -d "${TRANSCRIPTOME}" ]]; then
    echo "Error: transcriptome directory not found: ${TRANSCRIPTOME}" >&2
    exit 1
fi

if [[ ! -d "${FASTQ_DIR}" ]]; then
    echo "Error: FASTQ directory not found: ${FASTQ_DIR}" >&2
    exit 1
fi

if [[ ! -f "${IMAGE_FILE}" ]]; then
    echo "Error: histology image not found: ${IMAGE_FILE}" >&2
    exit 1
fi

if [[ "${SPACERANGER_BIN}" == */* ]]; then
    if [[ ! -x "${SPACERANGER_BIN}" ]]; then
        echo "Error: Space Ranger executable not found: ${SPACERANGER_BIN}" >&2
        exit 1
    fi
elif ! command -v "${SPACERANGER_BIN}" >/dev/null 2>&1; then
    echo "Error: '${SPACERANGER_BIN}' was not found in PATH." >&2
    echo "Set SPACERANGER_BIN to the Space Ranger executable if necessary." >&2
    exit 1
fi

if [[ ! "${LOCALCORES}" =~ ^[1-9][0-9]*$ ]]; then
    echo "Error: --localcores must be a positive integer." >&2
    exit 2
fi

if [[ ! "${LOCALMEM}" =~ ^[1-9][0-9]*$ ]]; then
    echo "Error: --localmem must be a positive integer." >&2
    exit 2
fi

mkdir -p "${OUTPUT_DIR}"
cd "${OUTPUT_DIR}"

echo "Running Space Ranger count"
echo "  run ID:        ${RUN_ID}"
echo "  sample:        ${SAMPLE_NAME}"
echo "  transcriptome: ${TRANSCRIPTOME}"
echo "  FASTQs:        ${FASTQ_DIR}"
echo "  image:         ${IMAGE_FILE}"
echo "  slide/area:    ${SLIDE_ID}/${CAPTURE_AREA}"
echo "  output:        ${OUTPUT_DIR}/${RUN_ID}"

"${SPACERANGER_BIN}" count \
    --id="${RUN_ID}" \
    --transcriptome="${TRANSCRIPTOME}" \
    --fastqs="${FASTQ_DIR}" \
    --sample="${SAMPLE_NAME}" \
    --image="${IMAGE_FILE}" \
    --slide="${SLIDE_ID}" \
    --area="${CAPTURE_AREA}" \
    --localcores="${LOCALCORES}" \
    --localmem="${LOCALMEM}"
