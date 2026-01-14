/*─────────────────────────────────────────────────────────────
  File: include/mesh_io.hpp
  Structured-zone utilities for FlowState BC_define
─────────────────────────────────────────────────────────────*/
#pragma once

#include "common.hpp"
#include <vector>
#include <string>

extern "C" {
    #include <cgnslib.h>          // CGNS C API
}

namespace fs {

/*-------------------------------------------------------------
  Zone – light wrapper around one CGNS Zone_t
-------------------------------------------------------------*/
struct Zone
{
    int         idx      = 0;           ///< 1-based zone index in CGNS file
    std::string name;                   ///< Zone_t name
    std::array<long long,3> vtxSize{};  ///< Vertex dimensions (Ni, Nj, Nk)

    /// Convenience: number of vertices along each axis
    long long ni() const { return vtxSize[0]; }
    long long nj() const { return vtxSize[1]; }
    long long nk() const { return vtxSize[2]; }

    /// Return PointRange on a whole-plane face
    PointRange face_range(FaceDir dir) const;
};

/*-------------------------------------------------------------
  Mesh – opens a CGNS file and collects structured Zone info
-------------------------------------------------------------*/
class Mesh
{
public:
    Mesh()  = default;
    ~Mesh();                           ///< ensures cg_close

    /// Open file read-write (modify) or read-only
    void open(const std::string& path, bool modify = true);
    void close();

    std::vector<Zone>&       zones()       { return zones_; }
    const std::vector<Zone>& zones() const { return zones_; }

    int file_id() const { return cgfile_; }
    int base_id() const { return cgbas_; }
    int cell_dim() const { return cell_dim_; }

private:
    int              cgfile_  = -1;
    int              cgbas_   =  1;        ///< assume first Base_t
    bool             isOpen_  = false;
    int              cell_dim_ = 0;
    std::vector<Zone> zones_;
};

} // namespace fs
