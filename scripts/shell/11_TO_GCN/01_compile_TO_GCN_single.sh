#!/usr/bin/env bash

# Download and compile the original single-time-series TO-GCN programs.
# No local absolute paths are used. The exact checked-out commit is recorded.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/../../.." && pwd)"
SOURCE_DIR="${PROJECT_ROOT}/external/TO-GCN"
BIN_DIR="${PROJECT_ROOT}/bin/TO-GCN"
LOG_DIR="${PROJECT_ROOT}/results/logs/11_TO_GCN"

native_path() {
    if command -v cygpath >/dev/null 2>&1; then
        cygpath -w "$1"
    else
        printf '%s\n' "$1"
    fi
}

SOURCE_DIR_NATIVE="$(native_path "${SOURCE_DIR}")"

TOGCN_REPOSITORY="${TOGCN_REPOSITORY:-https://github.com/petitmingchang/TO-GCN.git}"
TOGCN_REF="${TOGCN_REF:-master}"
CXX="${CXX:-g++}"

if ! command -v "${CXX}" >/dev/null 2>&1; then
    for rtools_compiler in \
        C:/rtools45/x86_64-w64-mingw32.static.posix/bin/g++.exe \
        C:/rtools44/x86_64-w64-mingw32.static.posix/bin/g++.exe \
        C:/rtools43/x86_64-w64-mingw32.static.posix/bin/g++.exe; do
        if [[ -x "${rtools_compiler}" ]]; then
            CXX="${rtools_compiler}"
            break
        fi
    done
fi

if ! command -v "${CXX}" >/dev/null 2>&1 && [[ ! -x "${CXX}" ]]; then
    echo "No C++ compiler was found. Install g++ or Rtools, or set CXX." >&2
    exit 1
fi

# An absolute Rtools g++ path does not automatically place its matching
# assembler and linker ahead of Cygwin utilities. Prepending the compiler
# directory prevents incompatible Cygwin binutils from being selected.
CXX_DIR="$(dirname "${CXX}")"
if [[ "${CXX_DIR}" != "." ]]; then
    if command -v cygpath >/dev/null 2>&1; then
        CXX_DIR_FOR_PATH="$(cygpath -u "${CXX_DIR}")"
    else
        CXX_DIR_FOR_PATH="${CXX_DIR}"
    fi
    export PATH="${CXX_DIR_FOR_PATH}:${PATH}"
fi

mkdir -p "${PROJECT_ROOT}/external" "${BIN_DIR}" "${LOG_DIR}"

if [[ ! -d "${SOURCE_DIR}/.git" ]]; then
    git clone "${TOGCN_REPOSITORY}" "${SOURCE_DIR_NATIVE}"
fi

if [[ "${TOGCN_REF}" != "master" ]]; then
    git -C "${SOURCE_DIR_NATIVE}" fetch --depth 1 origin "${TOGCN_REF}"
    git -C "${SOURCE_DIR_NATIVE}" checkout --detach FETCH_HEAD
fi

case "$(uname -s)" in
    MINGW*|MSYS*|CYGWIN*) executable_suffix=".exe" ;;
    *) executable_suffix="" ;;
esac

SOURCE_CODE_DIR="${SOURCE_DIR}/Single_Time-series_data/source_code"

for source_name in Cutoff_single GCN_single TO-GCN_single; do
    source_file="${SOURCE_CODE_DIR}/${source_name}.cpp"
    output_file="${BIN_DIR}/${source_name}${executable_suffix}"
    if [[ ! -f "${source_file}" ]]; then
        echo "Missing TO-GCN source file: ${source_file}" >&2
        exit 1
    fi
    "${CXX}" -O2 -std=gnu++11 \
        "$(native_path "${source_file}")" \
        -o "$(native_path "${output_file}")"
done

git -C "${SOURCE_DIR_NATIVE}" rev-parse HEAD > "${LOG_DIR}/TO_GCN_source_commit.txt"
"${CXX}" --version > "${LOG_DIR}/compiler_version.txt"

echo "Compiled TO-GCN programs in: ${BIN_DIR}"
echo "Source commit: $(cat "${LOG_DIR}/TO_GCN_source_commit.txt")"
