# Structured Mesh Connectivity & Boundary Conditions

bc_define is a preprocessing utility for structured CFD meshes.
It performs two primary tasks:

	1.	Detects and records 1-to-1 block connectivity
	2.	Applies user-defined boundary conditions to block faces or patches

The tool supports both CGNS (.cgns) and Plot3D (.x) multiblock structured mesh formats and is designed to be an early step in a CFD solver pipeline (e.g., FlowState).

# What bc_define Does

Given a structured mesh:

	•	Identifies abutting block faces that are geometrically coincident
	•	Verifies they are 1-to-1 conforming
	•	Writes connectivity metadata:
	  •	CGNS: GridConnectivity_t + GridConnectivityProperty_t
	  •	Plot3D: companion connectivity metadata file
	•	Applies boundary condition specifications from a user input file
	•	Writes all results back into the mesh file (CGNS) or into generated sidecar files (Plot3D)

## Requirements
- CMake 3.20+ or Make
- C++20 compiler
- CGNS (with HDF5). CMake can auto-fetch pinned CGNS + HDF5 when
  `AUTO_DEPS=ON`.

## Build
From the top-level directory:

### MacOS

macOS (CMake):
```
cmake -S . -B build
cmake --build build
```

macOS (Makefile):
```
make -C BC_define -f ../macos/BC_define.Makefile
```

### Linux

Linux (CMake):
```
sudo apt-get update
sudo apt-get install -y libhdf5-dev libcgns-dev
rm -rf build
cmake -S . -B build -DAUTO_DEPS=OFF
cmake --build build
```

Linux (Makefile):
```
make -C BC_define -f ../linux/BC_define.Makefile CGNS_INC=/path/include CGNS_LIB=/path/lib
```

Binaries are written to `BC_define/bin/`.

Note: CMake will try to find CGNS first. If it is not found and
`AUTO_DEPS=ON` (default), CMake downloads pinned CGNS + HDF5 into
the build tree. Set `-DAUTO_DEPS=OFF` and provide CGNS via
`CGNS_INCLUDE_DIR`/`CGNS_LIBRARY` or a `cgns` CMake package to use
a system install. AUTO_DEPS applies to CMake only; the Makefile
expects CGNS/HDF5 to be installed.

## Quick start
CGNS:
```
BC_define/bin/bc_define mesh.cgns bcdef.input --overwrite [--autowall | --autofarfield] --threads 16

BC_define/bin/bc_dump mesh.cgns
```

Plot3D:
```
BC_define/bin/bc_define mesh.x bcdef.input [--autowall | --autofarfield] --threads 16
```
Running bc_define with a Plot3D formatted mesh writes two files in the current working directory: `1to1s` and `bcs`. bc_dump only works with a .cgns file.

## bcdef.input format
Each non-comment line defines a BC family, a list of zones, and a face selector:
```
# comments start with '#'
# <bc_name>   <block_selection>   <location>
inflow   1,2   i-
outflow  3-4   i+
zerograd 2,4   j+
```
- Zones are 1-based indices.
- Face selector: `i-`, `i+`, `j-`, `j+`, `k-`, `k+`.
- Names are case-insensitive for type mapping, but the first spelling is kept
  for the family name (e.g., `Inflow` vs `inflow`).

## Boundary Condition Types

Supported BC names are normalized internally and mapped to CGNS enums:

| User Token  |  CGNS BC |
| -- | -- |
| inflow  |  BCInflow |
| outflow  | BCOutflow |
| zerograd  |  BCExtrapolate |
| symmetry  |  BCSymmetryPlane |
| wall  |  BCWall |
| wallviscous |  BCWallViscous |
| wallinviscid |  BCWallInviscid |
| farfield |  BCFarfield |
| extrapolate |  BCExtrapolate |
| neumann |  BCNeumann |
| dirichlet |  BCDirichlet |
| general |  BCGeneral |
| userdefined |  BCTypeUserDefined |

## Outputs
### CGNS
- Writes `ZoneGridConnectivity_t` and `ZoneBC_t` into the CGNS file.
- Boundary conditions are stored as `FamilySpecified` with a family node.
- `bc_dump` prints a human-readable summary.

### Plot3D
`bc_define` writes `1to1s` and `bcs` in the current working directory.

`1to1s` columns:
```
recv_zone donor_zone
r_i1 r_j1 r_k1 r_i2 r_j2 r_k2
d_i1 d_j1 d_k1 d_i2 d_j2 d_k2
t1 t2 t3
```
All ranges are inclusive, 1-based indices of nodes.

`bcs` columns:
```
zone family i1 j1 k1 i2 j2 k2
```
All ranges are inclusive, 1-based indices of cell edges.


## Behavior and assumptions
- Structured blocks only. 2D blocks are represented by `nk == 1`.
- 1-to-1 face connectivity only (matching grids).
- 3D blocks must be right-handed. Meshes with a left-handed zone are rejected.
- `--autowall` and `--autofarfield` are mutually exclusive.
- Overlapping user-specified BCs on the same face are treated as an error.

## Flags
- `--overwrite`    Purge existing connectivity/BCs in CGNS mesh before writing new ones.
- `--autowall`     Apply wall BCs to unassigned boundaries.
- `--autofarfield` Apply farfield BCs to unassigned boundaries.
- `--threads N`    Thread count for face-pair search.
- `--bench-io`     CGNS-only IO benchmark mode.
- `--bench-iter N` Iterations for `--bench-io`.

## Tests
CGNS tests:
```
./test/runtests_cgns.sh
```
Plot3D tests:
```
./test/runtests_p3d.sh
```

Expected outputs:
- `macos/expected/tests_correct_output_cgns`
- `macos/expected/tests_correct_output_plot3d`
- `linux/expected/` (generate by running tests on Linux)
