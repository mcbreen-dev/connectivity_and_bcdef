/*
  File: include/mesh_io.hpp

  CGNS structured-zone mesh wrapper used by bc_define and connectivity code.

  This header defines:
    - fs::Zone: lightweight per-zone metadata for a structured CGNS Zone_t
    - fs::Mesh: RAII-style opener/closer for a CGNS file and collector of zones

  Usage:
    - bc_define_main.cpp:
        Mesh::open(...) to open CGNS in modify/read mode
        Mesh::zones() to enumerate zones and pass Mesh to connectivity detector
    - connectivity_core.cpp + connectivity_*:
        Mesh::file_id(), Mesh::base_id(), Mesh::zones() to read coordinates and build faces
    - mesh_utils.cpp:
        Mesh::file_id(), Mesh::base_id(), Mesh::zones() to delete nodes / scan families
    - bc_define_cgns.cpp:
        Mesh::cell_dim() and Mesh::zones() to compute ranges and grid locations for BC writes
    - bc_dump.cpp:
        Mesh::open(..., false), Mesh::zones(), Mesh::file_id(), Mesh::base_id() to dump data

  Sizes / types:
    - CGNS file and base identifiers are stored as int (typically 32-bit).
    - Zone indices are 1-based ints matching CGNS zone numbering.
    - Vertex sizes are stored as long long (typically 64-bit) to avoid overflow in
      Ni*Nj*Nk computations in the connectivity code. (see common.hpp::IJK).
    - CGNS uses cgsize_t for indices; conversions occur in mesh_io.cpp and writers.
*/
#pragma once

#include "common.hpp"
#include <vector>
#include <string>

extern "C" {
    #include <cgnslib.h>          // CGNS C API; provides cg_open/cg_close/cg_zone_read/...
}

namespace fs {

/*-------------------------------------------------------------
  Zone – light wrapper around one CGNS Zone_t
-------------------------------------------------------------*/
struct Zone
{
    /*
      idx:
        - 1-based zone index as used by CGNS API calls (cg_zone_read, cg_nbocos, ...)
        - stored as int (typically 4 bytes)

      name:
        - Zone_t name read by cg_zone_read into a fixed-size char buffer
        - stored as std::string (owns a dynamically sized buffer)

      vtxSize:
        - vertex dimensions (Ni, Nj, Nk) for the structured zone
        - read from cg_zone_read's size[0][0..2] (cgsize_t)
        - cast to std::array<long long, 3> (24 bytes on LP64: 3 * 8)
        - for 2D meshes (cell_dim==2), Mesh::open forces Nk=1
    */
    int         idx      = 0;           ///< 1-based zone index in CGNS file
    std::string name;                   ///< Zone_t name
    std::array<long long,3> vtxSize{};  ///< Vertex dimensions (Ni, Nj, Nk)

    /// Convenience: number of vertices along each axis (long long, 64-bit)
    long long ni() const { return vtxSize[0]; }
    long long nj() const { return vtxSize[1]; }
    long long nk() const { return vtxSize[2]; }

    /*
      face_range(dir):
        - Returns a PointRange on the *vertex grid* describing the full face plane.
        - Output is 1-based indices.
        - Implemented in mesh_io.cpp.

      Usage:
        - NOT CURRENTLY USED: the connectivity code computes face ranges via
          other helpers (fill_point_range, face_center_range).
    */
    PointRange face_range(FaceDir dir) const;
};

/*-------------------------------------------------------------
  Mesh – opens a CGNS file and collects structured Zone info
-------------------------------------------------------------*/
class Mesh
{
public:
    Mesh()  = default;

    /*
      Destructor:
        - Calls close(), which calls cg_close if the file is open.
        - Implemented in mesh_io.cpp.

      Type notes:
        - Mesh is non-copyable by default due to std::vector and file handle semantics.
        - The codebase uses Mesh as a stack object in main() functions.
    */
    ~Mesh();    ///< ensures cg_close

    /*
      open(path, modify):
        - Opens CGNS file and populates zones_ with Zone metadata.

      Parameters:
        path   : filesystem path to .cgns file
        modify : if true, open with CG_MODE_MODIFY; else CG_MODE_READ

      Effects (mesh_io.cpp):
        - closes any previously open file
        - cg_open to get cgfile_
        - cg_base (assumes base index 1) and stores cell_dim_
        - cg_nzones to count zones
        - for each zone:
            cg_zone_read into name + size array
            stores Zone.idx, Zone.name, Zone.vtxSize
          If cell_dim_ == 2, forces Zone.vtxSize[2] = 1.

      Errors:
        - Uses CG_CALL to throw std::runtime_error on CGNS failures
    */
    void open(const std::string& path, bool modify = true);

    /*
      close():
        - If open, calls cg_close(cgfile_) and marks isOpen_=false
        - Safe to call multiple times
    */
    void close();

    /*
      zones():
        - Accessors for loaded zones.
        - Returned vector elements are stable unless zones_ is modified.
        - Used by virtually all subsystems (connectivity, BC writing, dumps).

      Returns:
        - reference to the internal vector<Zone>
        - zones are stored in CGNS zone order; zone idx is 1-based and equals (position+1)
    */
    std::vector<Zone>&       zones()       { return zones_; }
    const std::vector<Zone>& zones() const { return zones_; }

    /*
      file_id():
        - Returns cgfile_ (CGNS file handle).
        - Type: int (CGNS API uses int for file/base/zone handles).

      base_id():
        - Returns cgbas_ (Base_t index). The implementation assumes the first base (1).

      cell_dim():
        - Returns cell_dim_ read from cg_base_read. (2 or 3)
        - Used in bc_define_cgns.cpp to choose:
            * whether Nk is forced to 1 during open
            * BC GridLocation (EdgeCenter vs IFaceCenter/JFaceCenter/KFaceCenter)
            * 2D PointRange packing when writing/reading BCs and 1to1s
    */
    int file_id() const { return cgfile_; }
    int base_id() const { return cgbas_; }
    int cell_dim() const { return cell_dim_; }

private:
    /*
      cgfile_:
        - CGNS file handle (int). Initialized to -1.
        - Valid after cg_open succeeds.

      cgbas_:
        - Base index within CGNS file (int).
        - Set to 1 and not changed in the current codebase.

      isOpen_:
        - Tracks whether cgfile_ is currently open.

      cell_dim_:
        - Cell dimension read from Base_t: 2 or 3 in expected use cases.

      zones_:
        - Loaded zone metadata (Zone structs) collected during open().
    */
    int              cgfile_  = -1;
    int              cgbas_   =  1;        ///< assume first Base_t
    bool             isOpen_  = false;
    int              cell_dim_ = 0;
    std::vector<Zone> zones_;
};

} // namespace fs
