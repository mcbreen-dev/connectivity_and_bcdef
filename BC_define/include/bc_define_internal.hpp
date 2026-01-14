#pragma once

#include "common.hpp"
#include "mesh_io.hpp"
#include "plot3d_io.hpp"
#include "connectivity.hpp"

#include <algorithm>
#include <array>
#include <string>
#include <unordered_map>
#include <vector>

#include <cgnslib.h>

#ifndef CGNS_MAX_NAME_LENGTH
#define CGNS_MAX_NAME_LENGTH 32
#endif

namespace fs {
class Logger;
}

namespace fs::bc {

struct ZoneFaceKey {
    int zone = 0;
    fs::FaceDir face = fs::FaceDir::IMIN;

    bool operator==(const ZoneFaceKey& other) const
    {
        return zone == other.zone && face == other.face;
    }
};

struct ZoneFaceKeyHash {
    size_t operator()(const ZoneFaceKey& key) const
    {
        return (static_cast<size_t>(key.zone) << 3) ^
               static_cast<size_t>(key.face);
    }
};

struct BCSpec {
    std::string family;
    CGNS_ENUMT(BCType_t) type = CGNS_ENUMV(BCTypeUserDefined);
};

struct Plot3DBcEntry {
    int zone = 0;
    fs::FaceDir face = fs::FaceDir::IMIN;
    std::string family;
    std::array<long long, 6> range{};
};

std::string trim(const std::string& s);
std::string normalize_token(const std::string& s);
std::vector<int> parse_zone_list(const std::string& token);
CGNS_ENUMT(BCType_t) parse_bc_type(const std::string& token);
bool is_plot3d_mesh(const std::string& path);

fs::PointRange face_center_range(const fs::Zone& zone,
                                 fs::FaceDir face,
                                 const fs::IJK& vtxBegin,
                                 const fs::IJK& vtxEnd);
fs::PointRange face_center_range(const fs::Plot3DZone& zone,
                                 fs::FaceDir face,
                                 const fs::IJK& vtxBegin,
                                 const fs::IJK& vtxEnd);

CGNS_ENUMT(GridLocation_t) face_grid_location(fs::FaceDir face, int cell_dim);

template <typename T>
inline bool ranges_equal(const std::array<T, 6>& a, const std::array<T, 6>& b)
{
    return a == b;
}

template <typename T>
inline bool range_overlap(const std::array<T, 6>& a, const std::array<T, 6>& b,
                          std::array<T, 6>& out)
{
    for (int dim = 0; dim < 3; ++dim) {
        T lo = std::max(a[dim], b[dim]);
        T hi = std::min(a[dim + 3], b[dim + 3]);
        if (lo > hi)
            return false;
        out[dim] = lo;
        out[dim + 3] = hi;
    }
    return true;
}

template <typename T>
inline std::string range_to_string(const std::array<T, 6>& r)
{
    return "[" + std::to_string(r[0]) + "," + std::to_string(r[1]) + "," +
           std::to_string(r[2]) + "]-[" + std::to_string(r[3]) + "," +
           std::to_string(r[4]) + "," + std::to_string(r[5]) + "]";
}

std::string bc_patch_name(const std::string& family, fs::FaceDir face, int index);
bool range_face(const fs::Zone& zone,
                const std::array<cgsize_t, 6>& r,
                fs::FaceDir& out);
int ensure_family(int fn, int B,
                  const std::string& family,
                  CGNS_ENUMT(BCType_t) bc_type,
                  std::unordered_map<std::string, int>& cache);

std::unordered_map<ZoneFaceKey, BCSpec, ZoneFaceKeyHash>
parse_bcdef(const std::string& path);

void write_boundary_conditions(
    fs::Mesh& mesh,
    const std::vector<fs::BoundaryPatch>& patches,
    const std::unordered_map<ZoneFaceKey, BCSpec, ZoneFaceKeyHash>& specs,
    bool autowall,
    bool autofarfield);

void run_io_benchmark(fs::Mesh& mesh, fs::Logger& log, int iters);

void write_plot3d_1to1s(const std::string& path,
                        const std::vector<fs::Plot3DZone>& zones,
                        const std::vector<fs::ConnPatch>& conns);

std::vector<Plot3DBcEntry> build_plot3d_bcs(
    const std::vector<fs::Plot3DZone>& zones,
    const std::vector<fs::BoundaryPatch>& patches,
    const std::unordered_map<ZoneFaceKey, BCSpec, ZoneFaceKeyHash>& specs,
    bool autowall,
    bool autofarfield);

void write_plot3d_bcs(const std::string& path,
                      const std::vector<fs::Plot3DZone>& zones,
                      const std::vector<Plot3DBcEntry>& bcs);

} // namespace fs::bc
