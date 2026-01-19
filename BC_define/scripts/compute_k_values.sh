#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
ROOT_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
TOP_DIR="$(cd "${ROOT_DIR}/.." && pwd)"

OS_NAME="unknown"
case "$(uname -s)" in
    Darwin) OS_NAME="macos" ;;
    Linux) OS_NAME="linux" ;;
esac

MESH_SCRIPTS_DIR="${TOP_DIR}/${OS_NAME}/io_bench"

GEN_SCRIPT="${MESH_SCRIPTS_DIR}/generate_io_meshes.sh"
BENCH_SCRIPT="${MESH_SCRIPTS_DIR}/bench_io_k.sh"

if [[ ! -x "$GEN_SCRIPT" || ! -x "$BENCH_SCRIPT" ]]; then
    echo "Missing io_bench scripts for ${OS_NAME}." >&2
    echo "Expected: ${GEN_SCRIPT} and ${BENCH_SCRIPT}" >&2
    exit 1
fi

TMP_DIR="$(mktemp -d)"
cleanup() {
    rm -rf "$TMP_DIR"
}
trap cleanup EXIT

OUT_DIR="${TMP_DIR}" "$GEN_SCRIPT" >/dev/null
BENCH_OUTPUT="$(MESH_DIR="${TMP_DIR}" "$BENCH_SCRIPT")"

k2d="$(echo "$BENCH_OUTPUT" | awk -F': ' '/Suggested k_2d/ {print $2; exit}')"
k3d="$(echo "$BENCH_OUTPUT" | awk -F': ' '/Suggested k_3d/ {print $2; exit}')"

if [[ -z "${k2d}" || -z "${k3d}" ]]; then
    echo "Failed to extract k values" >&2
    exit 1
fi

cat <<EOF
#pragma once

namespace bcdef {
constexpr double K2D = ${k2d};
constexpr double K3D = ${k3d};
constexpr unsigned long long MEM_CAP_BYTES = 1ULL << 30;
}
EOF
