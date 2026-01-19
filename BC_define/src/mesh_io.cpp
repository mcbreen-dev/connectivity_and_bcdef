/*─────────────────────────────────────────────────────────────
  File: src/mesh_io.cpp

  CGNS mesh I/O and basic structured-zone helpers.

  This file implements:
    - Mesh::open / Mesh::close
        * Open a CGNS file in read-only or modify mode
        * Discover structured zones in the first base
        * Record per-zone vertex dimensions (Ni,Nj,Nk)
        * Normalize 2-D structured meshes so Nk == 1

  Scope and assumptions:
    - Only structured grids are handled here.
    - Only the first CGNS base (base_id() == 1) is used.
    - Zones are indexed 1-based externally (CGNS convention) but
      stored internally in a 0-based vector.

  Higher-level connectivity, BC creation, and overwrite logic live
  elsewhere (connectivity_* and mesh_utils.cpp).
─────────────────────────────────────────────────────────────*/
#include "mesh_io.hpp"
#include "common.hpp"

#ifndef CGNS_MAX_NAME_LENGTH
#define CGNS_MAX_NAME_LENGTH 32
#endif

namespace fs {

/*=====================================================================
  Mesh::open

  Open a CGNS file and populate the Mesh with structured zones.

  Parameters:
    path   : filesystem path to CGNS file
    modify : if true, open in CG_MODE_MODIFY (read/write);
             otherwise open read-only

  Behavior:
    - Closes any previously opened CGNS file
    - Opens the file with CGNS API
    - Reads the number of zones in the first base
    - Reads base cell dimension to determine 2-D vs 3-D
    - Reads each Zone_t name and vertex dimensions
    - Forces Nk = 1 for 2-D structured meshes

  Errors:
    - Any CGNS API failure throws via CG_CALL with a descriptive message
=====================================================================*/
void Mesh::open(const std::string& path, bool modify)
{
    if (isOpen_) close();

    int mode = modify ? CG_MODE_MODIFY : CG_MODE_READ;
    CG_CALL(cg_open(path.c_str(), mode, &cgfile_),
        "CGNS: cannot open file " + path);

    /*---------------------------------------------------------
      Collect structured zones from the first base
    ---------------------------------------------------------*/
    int nZones = 0;
    CG_CALL(cg_nzones(cgfile_, cgbas_, &nZones),
        "CGNS: cg_nzones failed");

    zones_.clear();
    zones_.reserve(static_cast<size_t>(nZones));

    /*---------------------------------------------------------
      Determine whether the mesh is 2-D or 3-D from the base
      cell dimension (cell_dim == 2 → structured 2-D mesh)
    ---------------------------------------------------------*/
    int celldim = 0, physdim = 0;
    char basename[CGNS_MAX_NAME_LENGTH + 1];
    CG_CALL(cg_base_read(cgfile_, cgbas_, basename, &celldim, &physdim),
        "CGNS: cg_base_read failed");
    cell_dim_ = celldim;
    bool is2D = (cell_dim_ == 2);


    /*---------------------------------------------------------
      Read each zone's name and vertex dimensions

      cg_zone_read returns a 3x3 array:
        size[0] = vertex dimensions
        size[1] = cell dimensions
        size[2] = boundary vertex dimensions

      Only size[0] (vertex dimensions) are used here.
    ---------------------------------------------------------*/
    for (int iz = 1; iz <= nZones; ++iz)
    {
        char zname[CGNS_MAX_NAME_LENGTH + 1] = ""; // 32 chars + null
        cgsize_t size[3][3];                          // 3x3 array

        /* Pass pointer to first element to satisfy CGNS API */
        CG_CALL(cg_zone_read(cgfile_, cgbas_, iz, zname, &size[0][0]),
            "CGNS: cg_zone_read failed");

        Zone z;
        z.idx  = iz;                  // 1-based CGNS zone index
        z.name = zname;
        z.vtxSize = { static_cast<long long>(size[0][0]),  // Ni
                      static_cast<long long>(size[0][1]),  // Nj
                      static_cast<long long>(size[0][2]) };// Nk

        /*-----------------------------------------------------
          Normalize 2-D structured zones:
          CGNS may store Nk > 1 even for 2-D grids; this tool
          treats all 2-D meshes as Nk == 1 consistently.
        -----------------------------------------------------*/
        if (is2D) {
            z.vtxSize[2] = 1;            // force Nk = 1 for 2-D structured zones
        }

        zones_.push_back(z);
    }

    isOpen_ = true;
}

/*=====================================================================
  Mesh::close

  Close the currently open CGNS file (if any).

  Safe to call multiple times; only acts when a file is open.
=====================================================================*/
void Mesh::close()
{
    if (isOpen_)
        CG_CALL(cg_close(cgfile_), "CGNS: cg_close failed");
    isOpen_ = false;
}

/*=====================================================================
  Mesh destructor

  Ensures the CGNS file is closed when the Mesh object is destroyed.
=====================================================================*/
Mesh::~Mesh() { close(); }

} // namespace fs
