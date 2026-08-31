#!/usr/bin/env bash

# Run the original single-time-series TO-GCN programs on the five ordered
# embryonic-leaf domains prepared by the R preprocessing script.
#
# Usage:
#   bash scripts/shell/11_TO_GCN/02_run_TO_GCN_single.sh [positive_PCC_cutoff]
#
# The recovered manuscript/PPT setting is 0.95. Users should inspect the
# Cutoff_single output and adjust this value for new datasets.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/../../.." && pwd)"
INPUT_DIR="${PROJECT_ROOT}/results/tables/11_TO_GCN/input"
OUTPUT_DIR="${PROJECT_ROOT}/results/tables/11_TO_GCN/original_TO_GCN"
BIN_DIR="${PROJECT_ROOT}/bin/TO-GCN"
LOG_DIR="${PROJECT_ROOT}/results/logs/11_TO_GCN"

POSITIVE_PCC_CUTOFF="${1:-0.95}"
NUMBER_OF_FEATURES=5

case "$(uname -s)" in
    MINGW*|MSYS*|CYGWIN*) executable_suffix=".exe" ;;
    *) executable_suffix="" ;;
esac

CUTOFF_BIN="${BIN_DIR}/Cutoff_single${executable_suffix}"
GCN_BIN="${BIN_DIR}/GCN_single${executable_suffix}"
TOGCN_BIN="${BIN_DIR}/TO-GCN_single${executable_suffix}"

TF_INPUT="${INPUT_DIR}/TF_expression.tsv"
ALL_GENE_INPUT="${INPUT_DIR}/All_gene_expression.tsv"
SEED_INPUT="${INPUT_DIR}/seeds.txt"

# The original C++ source stores input filenames in fixed 100-character
# buffers. Run from the output directory and pass short relative paths to avoid
# overflowing those legacy buffers on deeply nested Windows project paths.
TF_INPUT_RELATIVE="../input/TF_expression.tsv"
ALL_GENE_INPUT_RELATIVE="../input/All_gene_expression.tsv"
SEED_INPUT_RELATIVE="../input/seeds.txt"

for required_file in \
    "${CUTOFF_BIN}" "${GCN_BIN}" "${TOGCN_BIN}" \
    "${TF_INPUT}" "${ALL_GENE_INPUT}" "${SEED_INPUT}"; do
    if [[ ! -f "${required_file}" ]]; then
        echo "Required file not found: ${required_file}" >&2
        exit 1
    fi
done

mkdir -p "${OUTPUT_DIR}" "${LOG_DIR}"
cd "${OUTPUT_DIR}"

"${CUTOFF_BIN}" \
    "${NUMBER_OF_FEATURES}" \
    "${TF_INPUT_RELATIVE}" \
    > "${LOG_DIR}/Cutoff_single.log" 2>&1

"${GCN_BIN}" \
    "${NUMBER_OF_FEATURES}" \
    "${TF_INPUT_RELATIVE}" \
    "${ALL_GENE_INPUT_RELATIVE}" \
    "${POSITIVE_PCC_CUTOFF}" \
    > "${LOG_DIR}/GCN_single.log" 2>&1

"${TOGCN_BIN}" \
    "${NUMBER_OF_FEATURES}" \
    "${TF_INPUT_RELATIVE}" \
    "${SEED_INPUT_RELATIVE}" \
    "${POSITIVE_PCC_CUTOFF}" \
    > "${LOG_DIR}/TO-GCN_single.log" 2>&1

cat > "${OUTPUT_DIR}/run_parameters.tsv" <<EOF
parameter	value
number_of_features	${NUMBER_OF_FEATURES}
ordered_domains	SAM;P1_P2;P3;P4;P5
positive_PCC_cutoff	${POSITIVE_PCC_CUTOFF}
TF_input	${TF_INPUT_RELATIVE}
all_gene_input	${ALL_GENE_INPUT_RELATIVE}
seed_input	${SEED_INPUT_RELATIVE}
EOF

for expected_output in PCC_histogram.tsv C1+.csv Node_level.csv Node_relation.csv; do
    if [[ ! -s "${OUTPUT_DIR}/${expected_output}" ]]; then
        echo "Expected TO-GCN output is missing or empty: ${expected_output}" >&2
        exit 1
    fi
done

echo "TO-GCN completed. Outputs: ${OUTPUT_DIR}"
echo "Positive PCC cutoff: ${POSITIVE_PCC_CUTOFF}"
