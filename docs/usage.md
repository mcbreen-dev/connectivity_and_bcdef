# bc_define usage

This project provides two CLI tools:

- `bc_define`: detect 1-to-1 connectivity and apply BCs.
- `bc_dump`: dump CGNS connectivity and BC summaries.

## bc_define

CGNS:
```
bc_define mesh.cgns bcdef.input [--overwrite] [--autowall | --autofarfield]
                          [--threads N]
```

Plot3D:
```
bc_define mesh.x bcdef.input [--autowall | --autofarfield] [--threads N]
```

CGNS I/O benchmark (no modification):
```
bc_define mesh.cgns [bcdef.input] --bench-io [--bench-iter N]
```

### Arguments

- `mesh.cgns` / `mesh.x`: structured mesh to process. The format is auto-detected
  by extension (`.cgns` or `.x`).
- `bcdef.input`: text file describing boundary conditions (see below). Required
  for normal runs; not required for `--bench-io`.

### Flags

- `--overwrite`: delete existing `ZoneBC_t` and `ZoneGridConnectivity_t` nodes
  before rewriting connectivity and BCs. Also prunes unused `Family_t` nodes
  after the write.
- `--autowall`: create a `wall` BC on any boundary patch that does not already
  have a BC. When `--overwrite` is not set, autowall only fills boundaries that
  still lack any BC after existing BCs are considered.
- `--autofarfield`: create a `farfield` BC on any boundary patch that does not
  already have a BC. This is mutually exclusive with `--autowall`.
- `--threads N`: thread count for face-pair verification. `0` uses hardware
  concurrency.
- `--bench-io`: CGNS-only IO benchmark mode. Reads coordinates and prints timing
  and a heuristic rule of thumb.
- `--bench-iter N`: iteration count for `--bench-io` (default 3).

### Behavior summary

- Opens the CGNS file in modify mode (unless `--bench-io`).
- Optionally purges existing BC/connectivity nodes if `--overwrite` is set.
- Runs the connectivity detector to discover 1-to-1 links and boundary patches.
- Parses `bcdef.input` and writes BCs for any matching boundary patches.
- With `--autowall` or `--autofarfield`, assigns those BCs to remaining
  unassigned boundaries.
- Without `--overwrite`, it skips writing BCs or 1-to-1 links that already exist
  with the same name and point-range.
- Overlapping user-specified BCs on the same face are treated as an error.

### Output files

CGNS:
- The mesh is modified in place.
- `bc_define.log` is written alongside the mesh file.

Plot3D:
- Writes `1to1s` and `bcs` to the current working directory.
- `bc_define.log` is written alongside the mesh file.

The Plot3D output format is described in `docs/plot3d_formats.md`.

## bc_dump

```
bc_dump mesh.cgns
```

Prints a summary of zone sizes, 1-to-1 connectivity, and BCs. For BCs, the
output label uses the BC family name when present.

## bcdef.input format

Each non-empty, non-comment line has the form:

```
<bc_type> <zone_list> <face>
```

- Comments start with `#` and are ignored.
- `zone_list` is a comma- or dash-separated list, like `1,3,5` or `1-4`.
- `face` is one of `i-`, `i+`, `j-`, `j+`, `k-`, `k+`.

Example:

```
# square
inflow  1,2 i-
outflow 3-4 i+
zerograd 2,4 j+
```

The `bc_type` token is normalized and mapped to a CGNS `BCType_t`. Common
aliases include `inflow`, `outflow`, `symmetry`, `wall`, `wallviscous`,
`wallinviscid`, `farfield`, `zerograd`/`zerogradient`, `dirichlet`,
`neumann`, and `general`. Unknown tokens map to `BCTypeUserDefined`.

Family name behavior:

- The family name is derived from the first token on the line (e.g., `inflow`).
- `bc_define` writes BCs as `FamilySpecified` and attaches the family name.
- BC node names are generated as `<family>_<face>_<count>` to keep them unique.

## Build dependency note

CMake attempts to find CGNS (with HDF5). If it is not found and `AUTO_DEPS=ON`,
CMake will download and build pinned CGNS + HDF5 into the build tree. Set
`-DAUTO_DEPS=OFF` and provide `CGNS_INCLUDE_DIR`/`CGNS_LIBRARY` or a `cgns`
CMake package to use a system install. AUTO_DEPS applies to CMake only;
the Makefile expects CGNS/HDF5 to be installed. Use the OS-specific
Makefile in `macos/BC_define.Makefile` or `linux/BC_define.Makefile`.

Linux (system packages + CMake):
```
sudo apt-get update
sudo apt-get install -y libhdf5-dev libcgns-dev
rm -rf build
cmake -S . -B build -DAUTO_DEPS=OFF
cmake --build build
```

## Tests

- `test/runtests_cgns.sh` runs the CGNS suite and writes `test/runtests.log`.
- `test/runtests_p3d.sh` runs the Plot3D suite and writes `test/runtests.log`.
- `test/run1test.sh <n>` runs a single test and writes `test/test<n>.log`.
