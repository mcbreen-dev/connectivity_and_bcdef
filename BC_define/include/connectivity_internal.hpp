#pragma once

#include "common.hpp"
#include "connectivity.hpp"
#include "mesh_io.hpp"
#include "plot3d_io.hpp"

#include <cstddef>
#include <cstdint>
#include <unordered_map>
#include <vector>

namespace fs::conn {

struct Point {
    double x = 0.0;
    double y = 0.0;
    double z = 0.0;
};

struct FacePlane {
    int face = -1;
    long long dim1 = 0;
    long long dim2 = 0;
    std::vector<double> x, y, z;
};

struct VolumeCoords {
    long long ni = 0;
    long long nj = 0;
    long long nk = 0;
    std::vector<double> x, y, z;
};

struct FaceGrid {
    int uLen = 0;
    int vLen = 0;
    std::vector<int> index;
    std::vector<unsigned char> consumed;
};

struct FaceData {
    std::vector<int> block;
    std::vector<int> face;
    std::vector<int> u;
    std::vector<int> v;
    std::vector<int> i;
    std::vector<int> j;
    std::vector<int> k;
    std::vector<double> x;
    std::vector<double> y;
    std::vector<double> z;

    void reserve(size_t count)
    {
        block.reserve(count);
        face.reserve(count);
        u.reserve(count);
        v.reserve(count);
        i.reserve(count);
        j.reserve(count);
        k.reserve(count);
        x.reserve(count);
        y.reserve(count);
        z.reserve(count);
    }

    size_t size() const { return x.size(); }
    bool empty() const { return x.empty(); }

    int push(int block_v, int face_v, int u_v, int v_v,
             int i_v, int j_v, int k_v, const Point& center)
    {
        int idx = static_cast<int>(x.size());
        block.push_back(block_v);
        face.push_back(face_v);
        u.push_back(u_v);
        v.push_back(v_v);
        i.push_back(i_v);
        j.push_back(j_v);
        k.push_back(k_v);
        x.push_back(center.x);
        y.push_back(center.y);
        z.push_back(center.z);
        return idx;
    }
};

struct FaceBounds {
    bool valid = false;
    Point min{};
    Point max{};
};

struct FaceSummary {
    int face_id = 0;
    int block = 0;
    int face = 0;
    int uLen = 0;
    int vLen = 0;
    bool normal_valid = false;
    Point normal{};
    Point centroid{};
    FaceBounds bounds{};
};

struct FacePair {
    int a = 0;
    int b = 0;
};

struct MatchPair {
    int a = 0;
    int b = 0;
};

struct UVMatch {
    int au = 0;
    int av = 0;
    int bu = 0;
    int bv = 0;
    int a_idx = -1;
    int b_idx = -1;
};

struct NeighborInfo {
    int block = -1;
    int face = -1;
    int i = -1;
    int j = -1;
    int k = -1;
};

struct BoundaryRecord {
    int zone = 0;
    FaceDir face = FaceDir::IMIN;
    bool full = false;
    IJK begin{};
    IJK end{};
};

struct CellKey {
    std::int64_t ix = 0;
    std::int64_t iy = 0;
    std::int64_t iz = 0;

    bool operator==(const CellKey& other) const
    {
        return ix == other.ix && iy == other.iy && iz == other.iz;
    }
};

struct CellKeyHash {
    size_t operator()(const CellKey& key) const;
};

struct FaceHash {
    std::unordered_map<CellKey, std::vector<int>, CellKeyHash> cells;
};

struct FaceBuildResult {
    std::vector<FaceGrid> grids;
    std::vector<FaceBounds> bounds;
    std::vector<size_t> cell_start;
    std::vector<int> cell_count;
    std::vector<int> face_cells;
    FaceData faces;
    std::vector<Point> face_sum;
    std::vector<int> face_counts;
    size_t faces_loaded = 0;
};

struct HashBuildResult {
    std::vector<FaceHash> face_hashes;
    std::vector<FacePair> pairs;
    double hash_cell = 1.0;
};

size_t plane_index(const FacePlane& plane, long long a, long long b);
size_t face_grid_index(const FaceGrid& grid, int u, int v);
size_t volume_index(const VolumeCoords& vol, long long i, long long j, long long k);
size_t plot3d_index(const Plot3DZone& zone, long long i, long long j, long long k);

bool is_left_handed(int fn, int B, const Zone& z);
bool is_left_handed_plot3d(const Plot3DZone& z);

double distance_sq(const Point& a, const Point& b);
double distance_sq(double ax, double ay, double az, double bx, double by, double bz);
void expand_bounds(FaceBounds& bounds, const Point& p);
bool bounds_overlap(const FaceBounds& a, const FaceBounds& b, double tol);
Point cross(const Point& a, const Point& b);
double norm(const Point& a);
Point normalize(const Point& a);
double dot(const Point& a, const Point& b);

FaceDir face_dir_from_id(int face);
bool is_i_face(int face);
bool is_j_face(int face);
bool is_k_face(int face);

void face_dims(const Zone& z, int face, int& uLen, int& vLen);
void face_dims(const Plot3DZone& z, int face, int& uLen, int& vLen);
void face_cell_indices(const Zone& z, int face, int u, int v, int& i, int& j, int& k);
void face_cell_indices(const Plot3DZone& z, int face, int u, int v, int& i, int& j, int& k);

Point face_center_from_plane(const FacePlane& plane, const Zone& z, int face, int u, int v);
Point face_center_from_volume(const VolumeCoords& vol, const Zone& z, int face, int u, int v);
Point face_center_from_plot3d(const Plot3DZone& zone, int face, int u, int v);
FacePlane read_face_plane(int fn, int B, const Zone& z, int face);
VolumeCoords read_volume_coords(int fn, int B, const Zone& z);

bool should_write(int recvZone, int donorZone, int recvFace, int donorFace);
int signed_axis(int axis, int delta);
void set_normal_transform(int transform[3], int recvFace, int donorFace);

void fill_point_range(const Zone& z, int face, int uStart, int vStart, int uEnd, int vEnd,
                      IJK& begin, IJK& end);
void fill_point_range(const Plot3DZone& z, int face, int uStart, int vStart, int uEnd, int vEnd,
                      IJK& begin, IJK& end);

void compute_donor_range(const int transform[3], const IJK& recvBegin, const IJK& recvEnd,
                         const IJK& donorBeginIn, IJK& donorBegin, IJK& donorEnd);

std::int64_t quantize_coord(double v, double cell);
std::int64_t quantize_coord_offset(double v, double origin, double cell);
CellKey cell_key(const Point& p, double cell);
CellKey cell_key(double x, double y, double z, double cell);

FaceBuildResult build_faces_cgns(Mesh& mesh, bool show_progress);
FaceBuildResult build_faces_plot3d(const std::vector<Plot3DZone>& zones);
std::vector<FaceSummary> build_face_summaries(const FaceBuildResult& build);

HashBuildResult build_hash_candidates(const FaceBuildResult& build,
                                      const std::vector<FaceSummary>& summaries,
                                      double tol);

std::vector<MatchPair> collect_match_pairs(const FaceBuildResult& build,
                                           const std::vector<FaceSummary>& summaries,
                                           const HashBuildResult& hash,
                                           double tol,
                                           int threads);

void resolve_matches_cgns(const Mesh& mesh, Logger& log,
                          FaceBuildResult& build,
                          const std::vector<MatchPair>& matches,
                          std::vector<BoundaryRecord>& boundary_records,
                          std::vector<ConnPatch>& conn_patches);

void resolve_matches_plot3d(const std::vector<Plot3DZone>& zones, Logger& log,
                            FaceBuildResult& build,
                            const std::vector<MatchPair>& matches,
                            std::vector<BoundaryRecord>& boundary_records,
                            std::vector<ConnPatch>& conn_patches);

} // namespace fs::conn
