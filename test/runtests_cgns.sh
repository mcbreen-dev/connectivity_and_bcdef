#!/bin/bash

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
ROOT_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
cd "$SCRIPT_DIR" || exit 1

# Default logging: stdout to runtests.log. Use -v to also echo to terminal.
LOG_FILE="runtests.log"
if [ "${1:-}" = "-v" ]; then
    exec > >(tee "$LOG_FILE")
    shift
else
    exec > "$LOG_FILE"
fi
echo "PWD: $(pwd)"

# Define the path to the binaries (override with BC_DEFINE_DIR)
BIN_ROOT="${BC_DEFINE_DIR:-${ROOT_DIR}/BC_define}"
BIN="${BIN_ROOT}/bin/bc_define"
DUMP="${BIN_ROOT}/bin/bc_dump"

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
    "16 ${ROOT_DIR}/meshes/cgns/cube/cube3d_512M.cgns ${ROOT_DIR}/meshes/cgns/cube/bcdef.input_cube4"
    "17 ${ROOT_DIR}/meshes/cgns/tiny/tiny_square_4z.cgns ${ROOT_DIR}/meshes/cgns/square/bcdef.input_square1"
)

# Loop over each test case
for test in "${tests[@]}"; do

    num=$(echo "$test" | awk '{print $1}')
    mesh_file=$(echo "$test" | awk '{print $2}')
    input_file=$(echo "$test" | awk '{print $3}')
    mesh_basename=$(basename "$mesh_file")

    echo "======================"
    echo "       Test $num"
    echo "======================"

    # Remove any existing copy in current directory
    echo "Removing old mesh: $mesh_basename"
    rm -f "$mesh_basename"

    # Copy fresh version into the test directory and make writable
    echo "Copying mesh to current directory: $mesh_file"
    cp -f "$mesh_file" "./$mesh_basename"
    chmod u+w "./$mesh_basename"

    # Run bc_define
    echo "Running: $BIN $mesh_basename $input_file"
    $BIN "./$mesh_basename" "$input_file" "--overwrite" "--autowall" "--threads" "16"

    # Run bc_dump for every test
    echo "Running: $DUMP $mesh_basename"
    $DUMP "./$mesh_basename"

    echo "-------------------------------------------"
done

# Clean up mesh files and individual logs
rm -rf *.cgns bc_dump.log bc_define.log
