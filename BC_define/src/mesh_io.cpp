/*─────────────────────────────────────────────────────────────
  File: src/mesh_io.cpp
─────────────────────────────────────────────────────────────*/
#include "mesh_io.hpp"
#include "logger.hpp"
#include "common.hpp"

#ifndef CGNS_MAX_NAME_LENGTH
#define CGNS_MAX_NAME_LENGTH 32
#endif

namespace fs {

/*=========================  MESH ============================*/
void Mesh::open(const std::string& path, bool modify)
{
    if (isOpen_) close();

    int mode = modify ? CG_MODE_MODIFY : CG_MODE_READ;
    CG_CALL(cg_open(path.c_str(), mode, &cgfile_),
        "CGNS: cannot open file " + path);

    /* collect structured zones in first base */
    int nZones = 0;
    CG_CALL(cg_nzones(cgfile_, cgbas_, &nZones),
        "CGNS: cg_nzones failed");

    zones_.clear();
    zones_.reserve(static_cast<size_t>(nZones));

    /* Determine if mesh is 2D or 3D*/
    int celldim = 0, physdim = 0;
    char basename[CGNS_MAX_NAME_LENGTH + 1];
    CG_CALL(cg_base_read(cgfile_, cgbas_, basename, &celldim, &physdim),
        "CGNS: cg_base_read failed");
    cell_dim_ = celldim;
    bool is2D = (cell_dim_ == 2);

    for (int iz = 1; iz <= nZones; ++iz)
    {
        char zname[CGNS_MAX_NAME_LENGTH + 1] = ""; // 32 chars + null
        cgsize_t size[3][3];                          // 3x3 array

        /* Pass pointer to first element to satisfy cg_zone_read prototype */
        CG_CALL(cg_zone_read(cgfile_, cgbas_, iz, zname, &size[0][0]),
            "CGNS: cg_zone_read failed");

        Zone z;
        z.idx  = iz;
        z.name = zname;
        z.vtxSize = { static_cast<long long>(size[0][0]),  // Ni
                      static_cast<long long>(size[0][1]),  // Nj
                      static_cast<long long>(size[0][2]) };// Nk

        if (is2D) {
            z.vtxSize[2] = 1;            // force Nk = 1 for 2-D structured zones
        }

        zones_.push_back(z);
    }

    isOpen_ = true;
}

void Mesh::close()
{
    if (isOpen_)
        CG_CALL(cg_close(cgfile_), "CGNS: cg_close failed");
    isOpen_ = false;
}

Mesh::~Mesh() { close(); }

/*======================  ZONE HELPERS =======================*/
PointRange Zone::face_range(FaceDir dir) const
{
    PointRange pr;
    switch (dir)
    {
        case FaceDir::IMIN:
            pr.begin = {1, 1, 1};
            pr.end   = {1, nj(), nk()};
            break;

        case FaceDir::IMAX:
            pr.begin = {ni(), 1, 1};
            pr.end   = {ni(), nj(), nk()};
            break;

        case FaceDir::JMIN:
            pr.begin = {1, 1, 1};
            pr.end   = {ni(), 1, nk()};
            break;

        case FaceDir::JMAX:
            pr.begin = {1, nj(), 1};
            pr.end   = {ni(), nj(), nk()};
            break;

        case FaceDir::KMIN:
            pr.begin = {1, 1, 1};
            pr.end   = {ni(), nj(), 1};
            break;

        case FaceDir::KMAX:
            pr.begin = {1, 1, nk()};
            pr.end   = {ni(), nj(), nk()};
            break;
    }
    return pr;
}

} // namespace fs
