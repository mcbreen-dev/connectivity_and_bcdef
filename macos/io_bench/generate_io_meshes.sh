#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
TOP_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
SRC_DIR="${TOP_DIR}/BC_define/scripts/io_bench"
OUT_DIR="${OUT_DIR:-$SCRIPT_DIR}"
GEN_SRC="${SRC_DIR}/io_mesh_gen.cpp"
GEN_BIN="${SCRIPT_DIR}/io_mesh_gen"

CGNS_INC="${CGNS_INC:-/opt/homebrew/include}"
CGNS_LIB="${CGNS_LIB:-/opt/homebrew/lib}"
MAX_3D="${MAX_3D:-512}"

if [[ ! -x "$GEN_BIN" || "$GEN_SRC" -nt "$GEN_BIN" ]]; then
    clang++ -O2 -std=c++17 -I"$CGNS_INC" -L"$CGNS_LIB" \
        "$GEN_SRC" -lcgns -lhdf5 -o "$GEN_BIN"
fi

SPECS=$(cat <<'SPECS_EOF'
io2d_tiny_9        9    9    1   1 1 1
io2d_small_65      65   65   1   1 1 1
io2d_med_257       257  257  1   1 1 1
io2d_large_1025    1025 1025 1   1 1 1
io2d_slab_2049x33  2049 33   1   1 1 1
io2d_split_257_2x2 257  257  1   2 2 1
io2d_split_1025_4x4 1025 1025 1  4 4 1
io3d_tiny_9        9    9    9   1 1 1
io3d_small_33      33   33   33  1 1 1
io3d_med_129       129  129  129 1 1 1
io3d_large_257     257  257  257 1 1 1
io3d_max_512       512  512  512 1 1 1
io3d_slab_257x257x9 257 257  9   1 1 1
io3d_rod_512x17x17 512  17   17  1 1 1
io3d_split_129_2x2x2 129 129 129 2 2 2
io3d_split_257_4x4x4 257 257 257 4 4 4
SPECS_EOF
)

while read -r name ni nj nk sx sy sz; do
    [[ -z "${name}" ]] && continue
    if [[ "$nk" -gt 1 ]]; then
        if [[ "$ni" -gt "$MAX_3D" || "$nj" -gt "$MAX_3D" || "$nk" -gt "$MAX_3D" ]]; then
            echo "Skipping ${name} (${ni}x${nj}x${nk}) > MAX_3D=${MAX_3D}"
            continue
        fi
    fi

    out="${OUT_DIR}/${name}.cgns"
    echo "Generating ${out} (${ni}x${nj}x${nk}) split ${sx}x${sy}x${sz}"
    "$GEN_BIN" "$out" "$ni" "$nj" "$nk" "$sx" "$sy" "$sz"
done <<< "$SPECS"
