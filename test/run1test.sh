#!/bin/bash

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
ROOT_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
cd "$SCRIPT_DIR" || exit 1

if [ $# -ne 1 ]; then
    echo "Usage: $0 <test_number>"
    exit 1
fi

num="$1"

# Define the path to the binaries
BIN="${ROOT_DIR}/BC_define/bin/bc_define"
DUMP="${ROOT_DIR}/BC_define/bin/bc_dump"

# Define the list of tests: each entry is "mesh_file input_file"
tests=(
    "1 ${ROOT_DIR}/meshes/cgns/square/square2d_ij100.cgns ${ROOT_DIR}/meshes/cgns/square/bcdef.input_square1"
    "2 ${ROOT_DIR}/meshes/cgns/cube/cube3d_ijk100.cgns ${ROOT_DIR}/meshes/cgns/cube/bcdef.input_cube1"
    "3 ${ROOT_DIR}/meshes/cgns/square/square2d_z4ineg.cgns ${ROOT_DIR}/meshes/cgns/square/bcdef.input_square2"
    "4 ${ROOT_DIR}/meshes/cgns/square/square2d_z3ijneg.cgns ${ROOT_DIR}/meshes/cgns/square/bcdef.input_square3"
    "5 ${ROOT_DIR}/meshes/cgns/cube/cube3d_z6ikswap.cgns ${ROOT_DIR}/meshes/cgns/cube/bcdef.input_cube2"
    "6 ${ROOT_DIR}/meshes/cgns/cube/cube3d_z6lefthanded.cgns ${ROOT_DIR}/meshes/cgns/cube/bcdef.input_cube3"
    "7 ${ROOT_DIR}/meshes/cgns/rect/rect2d_3zones.cgns ${ROOT_DIR}/meshes/cgns/rect/bcdef.input_rect1"
    "8 ${ROOT_DIR}/meshes/cgns/rect/rect2d_4zones.cgns ${ROOT_DIR}/meshes/cgns/rect/bcdef.input_rect3"
    "9 ${ROOT_DIR}/meshes/cgns/rect/rect3d_4zones.cgns ${ROOT_DIR}/meshes/cgns/rect/bcdef.input_rect2"
    "10 ${ROOT_DIR}/meshes/cgns/rect/rect3d_6zones.cgns ${ROOT_DIR}/meshes/cgns/rect/bcdef.input_rect4"
    "11 ${ROOT_DIR}/meshes/cgns/skew/skew3d_6zones.cgns ${ROOT_DIR}/meshes/cgns/skew/bcdef.input_skew1"
    "12 ${ROOT_DIR}/meshes/cgns/blunt/blunt3d_1.cgns ${ROOT_DIR}/meshes/cgns/blunt/bcdef.input_blunt1"
    "13 ${ROOT_DIR}/meshes/cgns/blunt/blunt3d_2.cgns ${ROOT_DIR}/meshes/cgns/blunt/bcdef.input_blunt2"
    "14 ${ROOT_DIR}/meshes/cgns/omesh/omesh_1z.cgns ${ROOT_DIR}/meshes/cgns/omesh/bcdef.input_o"
    "15 ${ROOT_DIR}/meshes/cgns/omesh/dmesh_1z.cgns ${ROOT_DIR}/meshes/cgns/omesh/bcdef.input_b"
    "17 ${ROOT_DIR}/meshes/cgns/tiny/tiny_square_4z.cgns ${ROOT_DIR}/meshes/cgns/square/bcdef.input_square1"
)

test_line=""
for test in "${tests[@]}"; do
    test_num=$(echo "$test" | awk '{print $1}')
    if [ "$test_num" = "$num" ]; then
        test_line="$test"
        break
    fi
done

if [ -z "$test_line" ]; then
    echo "Unknown test number: $num"
    exit 1
fi

mesh_file=$(echo "$test_line" | awk '{print $2}')
input_file=$(echo "$test_line" | awk '{print $3}')
mesh_basename=$(basename "$mesh_file")

LOG_FILE="test${num}.log"
exec > "$LOG_FILE"
echo "PWD: $(pwd)"

# Remove any existing copy in current directory
rm -f "$mesh_basename"

# Copy fresh version into the test directory and make writable
cp -f "$mesh_file" "./$mesh_basename"
chmod u+w "./$mesh_basename"

# Run bc_define
$BIN "./$mesh_basename" "$input_file" "--autowall" "--overwrite"

# Run bc_dump
$DUMP "./$mesh_basename"

# Clean up logs
rm -rf bc_dump.log bc_define.log
