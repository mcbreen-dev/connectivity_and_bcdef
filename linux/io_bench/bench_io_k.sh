#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
TOP_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
SRC_DIR="${TOP_DIR}/BC_define/scripts/io_bench"
MESH_DIR="${MESH_DIR:-$SCRIPT_DIR}"

CXX="${CXX:-g++}"
CGNS_INC="${CGNS_INC:-/usr/include}"
CGNS_LIB="${CGNS_LIB:-/usr/lib}"

BC_DEFINE_BIN="${BC_DEFINE_BIN:-${SCRIPT_DIR}/io_bench}"
if [[ ! -x "$BC_DEFINE_BIN" || "${SRC_DIR}/io_bench.cpp" -nt "$BC_DEFINE_BIN" ]]; then
    "$CXX" -O2 -std=c++17 -I"$CGNS_INC" -L"$CGNS_LIB" \
        "${SRC_DIR}/io_bench.cpp" -lcgns -lhdf5 -o "$BC_DEFINE_BIN"
fi

BENCH_ITERS="${BENCH_ITERS:-3}"
TMP_DIR="$(mktemp -d)"
RESULTS="${TMP_DIR}/results.csv"

cleanup() {
    rm -rf "$TMP_DIR"
}
trap cleanup EXIT

echo "mesh,k_face,nk" > "$RESULTS"

shopt -s nullglob
for mesh in "${MESH_DIR}"/*.cgns; do
    base="$(basename "$mesh")"
    iters="$BENCH_ITERS"
    if [[ "$base" == *"max_512"* ]]; then
        iters=1
    fi

    echo "Benchmarking ${base} (iters=${iters})"
    out="$("${BC_DEFINE_BIN}" "$mesh" --bench-io --bench-iter "$iters")"

    k_face="$(echo "$out" | awk -F'k_face=' '/Overall:/ {split($2,a," "); print a[1]; exit}')"
    nk="$(echo "$out" | awk '/^Zone 1 / {
        n=split($NF,a,"x"); if (n==3) print a[3]; else print "";
        exit
    }')"

    if [[ -z "$k_face" || -z "$nk" ]]; then
        echo "Failed to parse ${base}" >&2
        continue
    fi

    echo "${base},${k_face},${nk}" >> "$RESULTS"
done

python3 - "$RESULTS" <<'PY'
import csv
import statistics
import sys

with open(sys.argv[1]) as f:
    rows = list(csv.DictReader(f))

def summarize(vals):
    vals = sorted(vals)
    if not vals:
        return None
    return {
        "min": min(vals),
        "median": statistics.median(vals),
        "max": max(vals),
        "count": len(vals),
    }

vals_2d = []
vals_3d = []
for r in rows:
    try:
        k = float(r["k_face"])
        nk = int(r["nk"])
    except (ValueError, TypeError):
        continue
    if nk == 1:
        vals_2d.append(k)
    else:
        vals_3d.append(k)

sum_2d = summarize(vals_2d)
sum_3d = summarize(vals_3d)

print("\nSummary (k_face):")
if sum_2d:
    print("  2D: count={count} min={min:.3f} median={median:.3f} max={max:.3f}".format(**sum_2d))
    print("  Suggested k_2d (median): {:.3f}".format(sum_2d["median"]))
else:
    print("  2D: no data")

if sum_3d:
    print("  3D: count={count} min={min:.3f} median={median:.3f} max={max:.3f}".format(**sum_3d))
    print("  Suggested k_3d (median): {:.3f}".format(sum_3d["median"]))
else:
    print("  3D: no data")
PY
