# Module map

This project is organized around the BC_define tool and its supporting utilities.

## Top level
- BC_define/     - C++ source, headers, build artifacts, and scripts
- build/         - Top-level CMake build output
- docs/          - User and developer documentation
- macos/         - macOS-specific Makefile, IO bench scripts, and expected outputs
- linux/         - Linux-specific Makefile, IO bench scripts, and expected outputs
- meshes/        - Test meshes and bcdef inputs
- test/          - Test scripts and working files

## Docs
- docs/usage.md
  - User-facing CLI usage and file formats
- docs/algorithm.md
  - Connectivity and BC-writing algorithm overview
- docs/cgns_conventions.md
  - CGNS conventions and metadata layout
- docs/plot3d_formats.md
  - Plot3D output file formats (`1to1s` and `bcs`)
- docs/testing.md
  - Test scripts and expected outputs
- docs/bc_dump.md
  - bc_dump output format
- docs/module_map.md
  - This module map

## BC_define
### Entry points
- BC_define/src/bc_define_main.cpp
  - CLI parsing, flag handling, high-level flow
- BC_define/src/bc_dump.cpp
  - CGNS metadata dumper for connectivity/BC inspection

### bc_define modules
- BC_define/src/bc_define_cgns.cpp
  - CGNS mesh I/O, BC writing, overwrite/autowall/autofarfield behavior
- BC_define/src/bc_define_plot3d.cpp
  - Plot3D output writers for 1to1s and bcs
- BC_define/src/bc_define_bcspec.cpp
  - bcdef input parsing + validation
- BC_define/src/bc_define_helpers.cpp
  - Small string/path/face helpers
- BC_define/include/bc_define_internal.hpp
  - Internal structs and shared helpers for bc_define

### Connectivity engine
- BC_define/src/connectivity_core.cpp
  - ConnectivityDetector::run/run_collect/run_plot3d
- BC_define/src/connectivity_face_build.cpp
  - Face extraction and face metadata (bounds/centers/normals)
- BC_define/src/connectivity_hash.cpp
  - Spatial hashing and candidate face pair generation
- BC_define/src/connectivity_verify.cpp
  - Pair verification + transform resolution
- BC_define/src/connectivity_utils.cpp
  - Shared math helpers and indexing utilities
- BC_define/include/connectivity_internal.hpp
  - Internal connectivity data structures and helper signatures

### CGNS writers
- BC_define/src/connectivity_writer.cpp
  - Writes 1-to-1 connectivity into CGNS
- BC_define/include/connectivity_writer.hpp
  - ConnPatch definition and write_1to1 declaration

### Mesh I/O
- BC_define/src/mesh_io.cpp / BC_define/include/mesh_io.hpp
  - CGNS mesh open/read helpers
- BC_define/src/plot3d_io.cpp / BC_define/include/plot3d_io.hpp
  - Plot3D .x reader
- BC_define/src/mesh_utils.cpp / BC_define/include/mesh_utils.hpp
  - CGNS cleanup helpers (overwrite, prune families)

### Logging
- BC_define/src/logger.cpp / BC_define/include/logger.hpp
  - Minimal logging wrapper writing bc_define.log or bc_dump.log

### Build and scripts
- BC_define/CMakeLists.txt
  - CMake build targets
- BC_define/Makefile
  - Stub Makefile that points to OS-specific Makefiles
- macos/BC_define.Makefile
  - macOS Makefile build targets
- linux/BC_define.Makefile
  - Linux Makefile build targets
- macos/cmake/Defaults.cmake
  - macOS CMake defaults (Homebrew CGNS/HDF5 paths)
- linux/cmake/Defaults.cmake
  - Linux CMake defaults (common CGNS/HDF5 paths)
- BC_define/scripts/compute_k_values.sh
  - Computes k2d/k3d used for IO strategy selection
- BC_define/scripts/io_bench/
  - IO bench source files used by OS-specific scripts
- macos/io_bench/
  - macOS IO bench scripts and binaries
- linux/io_bench/
  - Linux IO bench scripts and binaries

## Tests
- test/runtests_cgns.sh
  - CGNS integration tests
- test/runtests_p3d.sh
  - Plot3D integration tests
- test/run1test.sh
  - Runs a single test case
- macos/expected/tests_correct_output_cgns
  - Expected CGNS test output (macOS)
- macos/expected/tests_correct_output_plot3d
  - Expected Plot3D test output (macOS)
- linux/expected/
  - Place Linux expected outputs here
