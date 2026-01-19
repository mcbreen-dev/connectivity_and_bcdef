/*
  File: include/plot3d_io.hpp

  Plot3D (structured multi-block) mesh representation and reader declaration.

  This header defines:
    - fs::Plot3DZone: in-memory storage for one Plot3D block's vertex grid and coordinates
    - fs::read_plot3d(): reads a Plot3D file and returns a vector of Plot3DZone

  Usage across the codebase:
    - bc_define_main.cpp:
        - Detects Plot3D input by file extension ".x" (bc_define_helpers.cpp::is_plot3d_mesh)
        - Calls read_plot3d(meshPath) to load zones
        - Passes zones into ConnectivityDetector::run_plot3d(...)
        - Uses zones again when writing "1to1s" and "bcs" text outputs

    - connectivity_face_build.cpp + connectivity_utils.cpp:
        - Access Plot3DZone::ni/nj/nk and x/y/z storage to compute face centers and handedness

  Data layout / sizes:
    - vtxSize holds (Ni, Nj, Nk) as long long (3 * 8 bytes on LP64).
    - x/y/z store vertex coordinates as double (8 bytes each).
    - The coordinate vectors are sized to Ni*Nj*Nk.
    - plot3d_index(...) in connectivity_utils.cpp uses 1-based (i,j,k) inputs and computes
      a 0-based linear index into these vectors.

  Reader implementation summary (src/plot3d_io.cpp):
    - Treats file as Fortran unformatted records:
        [uint32 leading_bytes][payload...][uint32 trailing_bytes]
      and checks leading == trailing.
    - Reads:
        1) block_count as int32 scalar record
        2) a record containing block_count * 3 int32 dimensions
        3) for each block, a record containing (Ni*Nj*Nk*3) doubles
           laid out as all X, then all Y, then all Z.
*/
#pragma once

#include <array>
#include <cstdint>
#include <string>
#include <vector>

namespace fs {

struct Plot3DZone {
    /*
      idx:
        - 1-based block index assigned by read_plot3d (block+1)
        - type: int (typically 4 bytes)

      name:
        - generated name "block-N" in read_plot3d
        - type: std::string

      vtxSize:
        - vertex dimensions (Ni, Nj, Nk) stored as long long
        - values originate from Plot3D int32 dimensions and are widened to long long

      x/y/z:
        - coordinate arrays, each length = Ni*Nj*Nk
        - element type double (8 bytes)
        - indexing convention used by connectivity code:
            i in [1..Ni], j in [1..Nj], k in [1..Nk]
            linear index = (k-1)*Ni*Nj + (j-1)*Ni + (i-1)
    */
    int idx = 0; // 1-based
    std::string name;
    std::array<long long, 3> vtxSize{};
    std::vector<double> x;
    std::vector<double> y;
    std::vector<double> z;

    // Vertex counts along each axis (64-bit)
    long long ni() const { return vtxSize[0]; }
    long long nj() const { return vtxSize[1]; }
    long long nk() const { return vtxSize[2]; }
};

/*
  read_plot3d(mesh_path)

  Returns:
    - std::vector<Plot3DZone>, one entry per block in the file

  Throws:
    - std::runtime_error on file open failure, record marker mismatch,
      or unexpected record sizes.

  Called in:
    - bc_define_main.cpp (Plot3D workflow)
*/
std::vector<Plot3DZone> read_plot3d(const std::string& mesh_path);

} // namespace fs
