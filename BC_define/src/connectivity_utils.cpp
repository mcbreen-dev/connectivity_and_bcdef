//─────────────────────────────────────────────────────────────
// File: src/connectivity_utils.cpp
// Shared math and indexing utilities for connectivity detection
//─────────────────────────────────────────────────────────────
#include "connectivity_internal.hpp"

#include <algorithm>
#include <cmath>

extern "C" {
    #include <cgnslib.h>
}

namespace fs::conn {

size_t plane_index(const FacePlane& plane, long long a, long long b)
{
    return static_cast<size_t>((b - 1) * plane.dim1 + (a - 1));
}

size_t face_grid_index(const FaceGrid& grid, int u, int v)
{
    return static_cast<size_t>(v) * grid.uLen + u;
}

size_t volume_index(const VolumeCoords& vol, long long i, long long j, long long k)
{
    return static_cast<size_t>((k - 1) * vol.ni * vol.nj + (j - 1) * vol.ni + (i - 1));
}

size_t plot3d_index(const Plot3DZone& zone, long long i, long long j, long long k)
{
    return static_cast<size_t>((k - 1) * zone.ni() * zone.nj() + (j - 1) * zone.ni() + (i - 1));
}

bool is_left_handed(int fn, int B, const Zone& z)
{
    static constexpr const char* cname[3] =
        {"CoordinateX", "CoordinateY", "CoordinateZ"};

    auto node = [&](long long i, long long j, long long k)
    {
        double X, Y, Z;
        cgsize_t s[3] = {cgsize_t(i), cgsize_t(j), cgsize_t(k)};
        cgsize_t e[3] = {s[0], s[1], s[2]};
        CG_CALL(cg_coord_read(fn, B, z.idx, cname[0], CGNS_ENUMV(RealDouble), s, e, &X),
                "cg_coord_read failed");
        CG_CALL(cg_coord_read(fn, B, z.idx, cname[1], CGNS_ENUMV(RealDouble), s, e, &Y),
                "cg_coord_read failed");
        CG_CALL(cg_coord_read(fn, B, z.idx, cname[2], CGNS_ENUMV(RealDouble), s, e, &Z),
                "cg_coord_read failed");
        return Point{X, Y, Z};
    };

    Point p000 = node(1, 1, 1);
    Point di   = node(2, 1, 1);
    Point dj   = node(1, 2, 1);
    Point dk   = node(1, 1, (z.nk() > 1 ? 2 : 1));

    di.x -= p000.x; di.y -= p000.y; di.z -= p000.z;
    dj.x -= p000.x; dj.y -= p000.y; dj.z -= p000.z;
    dk.x -= p000.x; dk.y -= p000.y; dk.z -= p000.z;

    double det = di.x * (dj.y * dk.z - dj.z * dk.y)
               - di.y * (dj.x * dk.z - dj.z * dk.x)
               + di.z * (dj.x * dk.y - dj.y * dk.x);
    return det < 0.0;
}

bool is_left_handed_plot3d(const Plot3DZone& z)
{
    auto node = [&](long long i, long long j, long long k)
    {
        size_t idx = plot3d_index(z, i, j, k);
        return Point{z.x[idx], z.y[idx], z.z[idx]};
    };

    Point p000 = node(1, 1, 1);
    Point di   = node(2, 1, 1);
    Point dj   = node(1, 2, 1);
    Point dk   = node(1, 1, (z.nk() > 1 ? 2 : 1));

    di.x -= p000.x; di.y -= p000.y; di.z -= p000.z;
    dj.x -= p000.x; dj.y -= p000.y; dj.z -= p000.z;
    dk.x -= p000.x; dk.y -= p000.y; dk.z -= p000.z;

    double det = di.x * (dj.y * dk.z - dj.z * dk.y)
               - di.y * (dj.x * dk.z - dj.z * dk.x)
               + di.z * (dj.x * dk.y - dj.y * dk.x);
    return det < 0.0;
}

double distance_sq(const Point& a, const Point& b)
{
    double dx = a.x - b.x;
    double dy = a.y - b.y;
    double dz = a.z - b.z;
    return dx * dx + dy * dy + dz * dz;
}

double distance_sq(double ax, double ay, double az,
                   double bx, double by, double bz)
{
    double dx = ax - bx;
    double dy = ay - by;
    double dz = az - bz;
    return dx * dx + dy * dy + dz * dz;
}

void expand_bounds(FaceBounds& bounds, const Point& p)
{
    if (!bounds.valid) {
        bounds.valid = true;
        bounds.min = p;
        bounds.max = p;
        return;
    }
    bounds.min.x = std::min(bounds.min.x, p.x);
    bounds.min.y = std::min(bounds.min.y, p.y);
    bounds.min.z = std::min(bounds.min.z, p.z);
    bounds.max.x = std::max(bounds.max.x, p.x);
    bounds.max.y = std::max(bounds.max.y, p.y);
    bounds.max.z = std::max(bounds.max.z, p.z);
}

bool bounds_overlap(const FaceBounds& a, const FaceBounds& b, double tol)
{
    if (!a.valid || !b.valid)
        return false;
    if (a.min.x > b.max.x + tol || b.min.x > a.max.x + tol) return false;
    if (a.min.y > b.max.y + tol || b.min.y > a.max.y + tol) return false;
    if (a.min.z > b.max.z + tol || b.min.z > a.max.z + tol) return false;
    return true;
}

Point cross(const Point& a, const Point& b)
{
    return Point{
        a.y * b.z - a.z * b.y,
        a.z * b.x - a.x * b.z,
        a.x * b.y - a.y * b.x
    };
}

double norm(const Point& a)
{
    return std::sqrt(a.x * a.x + a.y * a.y + a.z * a.z);
}

Point normalize(const Point& a)
{
    double n = norm(a);
    if (n == 0.0)
        return Point{0, 0, 0};
    return Point{a.x / n, a.y / n, a.z / n};
}

double dot(const Point& a, const Point& b)
{
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

size_t CellKeyHash::operator()(const CellKey& key) const
{
    size_t h = static_cast<size_t>(key.ix);
    h ^= static_cast<size_t>(key.iy + 0x9e3779b97f4a7c15ULL) + (h << 6) + (h >> 2);
    h ^= static_cast<size_t>(key.iz + 0x9e3779b97f4a7c15ULL) + (h << 6) + (h >> 2);
    return h;
}

std::int64_t quantize_coord(double v, double cell)
{
    return static_cast<std::int64_t>(std::floor(v / cell));
}

std::int64_t quantize_coord_offset(double v, double origin, double cell)
{
    return static_cast<std::int64_t>(std::floor((v - origin) / cell));
}

CellKey cell_key(const Point& p, double cell)
{
    return CellKey{
        quantize_coord(p.x, cell),
        quantize_coord(p.y, cell),
        quantize_coord(p.z, cell)
    };
}

CellKey cell_key(double x, double y, double z, double cell)
{
    return CellKey{
        quantize_coord(x, cell),
        quantize_coord(y, cell),
        quantize_coord(z, cell)
    };
}

FaceDir face_dir_from_id(int face)
{
    switch (face) {
        case 0: return FaceDir::IMIN;
        case 1: return FaceDir::IMAX;
        case 2: return FaceDir::JMIN;
        case 3: return FaceDir::JMAX;
        case 4: return FaceDir::KMIN;
        case 5: return FaceDir::KMAX;
    }
    return FaceDir::IMIN;
}

bool is_i_face(int face) { return face == 0 || face == 1; }
bool is_j_face(int face) { return face == 2 || face == 3; }
bool is_k_face(int face) { return face == 4 || face == 5; }

void face_dims(const Zone& z, int face, int& uLen, int& vLen)
{
    if (z.nk() == 1) {
        if (face == 0 || face == 1) { // iMin/iMax
            uLen = static_cast<int>(z.nj() - 1);
            vLen = 1;
        } else if (face == 2 || face == 3) { // jMin/jMax
            uLen = static_cast<int>(z.ni() - 1);
            vLen = 1;
        } else {
            uLen = 0;
            vLen = 0;
        }
        return;
    }

    if (face == 0 || face == 1) {
        uLen = static_cast<int>(z.nj() - 1);
        vLen = static_cast<int>(z.nk() - 1);
    } else if (face == 2 || face == 3) {
        uLen = static_cast<int>(z.ni() - 1);
        vLen = static_cast<int>(z.nk() - 1);
    } else {
        uLen = static_cast<int>(z.ni() - 1);
        vLen = static_cast<int>(z.nj() - 1);
    }
}

void face_dims(const Plot3DZone& z, int face, int& uLen, int& vLen)
{
    if (z.nk() == 1) {
        if (face == 0 || face == 1) { // iMin/iMax
            uLen = static_cast<int>(z.nj() - 1);
            vLen = 1;
        } else if (face == 2 || face == 3) { // jMin/jMax
            uLen = static_cast<int>(z.ni() - 1);
            vLen = 1;
        } else {
            uLen = 0;
            vLen = 0;
        }
        return;
    }

    if (face == 0 || face == 1) {
        uLen = static_cast<int>(z.nj() - 1);
        vLen = static_cast<int>(z.nk() - 1);
    } else if (face == 2 || face == 3) {
        uLen = static_cast<int>(z.ni() - 1);
        vLen = static_cast<int>(z.nk() - 1);
    } else {
        uLen = static_cast<int>(z.ni() - 1);
        vLen = static_cast<int>(z.nj() - 1);
    }
}

void face_cell_indices(const Zone& z, int face, int u, int v, int& i, int& j, int& k)
{
    bool is2D = (z.nk() == 1);

    switch (face) {
        case 0: // iMin
            i = 1;
            j = u + 1;
            k = is2D ? 1 : (v + 1);
            break;
        case 1: // iMax
            i = static_cast<int>(z.ni());
            j = u + 1;
            k = is2D ? 1 : (v + 1);
            break;
        case 2: // jMin
            j = 1;
            i = u + 1;
            k = is2D ? 1 : (v + 1);
            break;
        case 3: // jMax
            j = static_cast<int>(z.nj());
            i = u + 1;
            k = is2D ? 1 : (v + 1);
            break;
        case 4: // kMin
            k = 1;
            i = u + 1;
            j = v + 1;
            break;
        case 5: // kMax
            k = static_cast<int>(z.nk());
            i = u + 1;
            j = v + 1;
            break;
    }
}

void face_cell_indices(const Plot3DZone& z, int face, int u, int v, int& i, int& j, int& k)
{
    bool is2D = (z.nk() == 1);

    switch (face) {
        case 0: // iMin
            i = 1;
            j = u + 1;
            k = is2D ? 1 : (v + 1);
            break;
        case 1: // iMax
            i = static_cast<int>(z.ni());
            j = u + 1;
            k = is2D ? 1 : (v + 1);
            break;
        case 2: // jMin
            j = 1;
            i = u + 1;
            k = is2D ? 1 : (v + 1);
            break;
        case 3: // jMax
            j = static_cast<int>(z.nj());
            i = u + 1;
            k = is2D ? 1 : (v + 1);
            break;
        case 4: // kMin
            k = 1;
            i = u + 1;
            j = v + 1;
            break;
        case 5: // kMax
            k = static_cast<int>(z.nk());
            i = u + 1;
            j = v + 1;
            break;
    }
}

bool should_write(int recvZone, int donorZone, int recvFace, int donorFace)
{
    if (recvZone < donorZone) return true;
    if (recvZone > donorZone) return false;
    return recvFace < donorFace;
}

int signed_axis(int axis, int delta)
{
    return (delta >= 0) ? axis : -axis;
}

void set_normal_transform(int transform[3], int recvFace, int donorFace)
{
    int recvId = recvFace + 1;
    int donorId = donorFace + 1;
    int normalIndexSign = ((recvId + donorId) % 2 == 0) ? -1 : 1;

    int donorAxis = 1;
    if (is_j_face(donorFace)) donorAxis = 2;
    else if (is_k_face(donorFace)) donorAxis = 3;

    if (is_i_face(recvFace)) {
        transform[0] = normalIndexSign * donorAxis;
    } else if (is_j_face(recvFace)) {
        transform[1] = normalIndexSign * donorAxis;
    } else {
        transform[2] = normalIndexSign * donorAxis;
    }
}

void fill_point_range(const Zone& z, int face, int uStart, int vStart, int uEnd, int vEnd,
                      IJK& begin, IJK& end)
{
    long long u0 = uStart + 1;
    long long u1 = uEnd + 2;
    long long v0 = vStart + 1;
    long long v1 = vEnd + 2;
    bool is2D = (z.nk() == 1);

    switch (face) {
        case 0: // iMin
            begin = {1, u0, is2D ? 1 : v0};
            end   = {1, u1, is2D ? 1 : v1};
            break;
        case 1: // iMax
            begin = {z.ni(), u0, is2D ? 1 : v0};
            end   = {z.ni(), u1, is2D ? 1 : v1};
            break;
        case 2: // jMin
            begin = {u0, 1, is2D ? 1 : v0};
            end   = {u1, 1, is2D ? 1 : v1};
            break;
        case 3: // jMax
            begin = {u0, z.nj(), is2D ? 1 : v0};
            end   = {u1, z.nj(), is2D ? 1 : v1};
            break;
        case 4: // kMin
            begin = {u0, v0, 1};
            end   = {u1, v1, 1};
            break;
        case 5: // kMax
            begin = {u0, v0, z.nk()};
            end   = {u1, v1, z.nk()};
            break;
    }
}

void fill_point_range(const Plot3DZone& z, int face, int uStart, int vStart, int uEnd, int vEnd,
                      IJK& begin, IJK& end)
{
    long long u0 = uStart + 1;
    long long u1 = uEnd + 2;
    long long v0 = vStart + 1;
    long long v1 = vEnd + 2;
    bool is2D = (z.nk() == 1);

    switch (face) {
        case 0: // iMin
            begin = {1, u0, is2D ? 1 : v0};
            end   = {1, u1, is2D ? 1 : v1};
            break;
        case 1: // iMax
            begin = {z.ni(), u0, is2D ? 1 : v0};
            end   = {z.ni(), u1, is2D ? 1 : v1};
            break;
        case 2: // jMin
            begin = {u0, 1, is2D ? 1 : v0};
            end   = {u1, 1, is2D ? 1 : v1};
            break;
        case 3: // jMax
            begin = {u0, z.nj(), is2D ? 1 : v0};
            end   = {u1, z.nj(), is2D ? 1 : v1};
            break;
        case 4: // kMin
            begin = {u0, v0, 1};
            end   = {u1, v1, 1};
            break;
        case 5: // kMax
            begin = {u0, v0, z.nk()};
            end   = {u1, v1, z.nk()};
            break;
    }
}

void compute_donor_range(const int transform[3], const IJK& recvBegin, const IJK& recvEnd,
                         const IJK& donorBeginIn, IJK& donorBegin, IJK& donorEnd)
{
    auto sgn = [](int x) { return x >= 0 ? 1 : -1; };
    auto delta = [](int x, int axis) { return (std::abs(x) == axis) ? 1 : 0; };

    int a = transform[0];
    int b = transform[1];
    int c = transform[2];

    int T[3][3];
    T[0][0] = sgn(a) * delta(a, 1);
    T[0][1] = sgn(b) * delta(b, 1);
    T[0][2] = sgn(c) * delta(c, 1);

    T[1][0] = sgn(a) * delta(a, 2);
    T[1][1] = sgn(b) * delta(b, 2);
    T[1][2] = sgn(c) * delta(c, 2);

    T[2][0] = sgn(a) * delta(a, 3);
    T[2][1] = sgn(b) * delta(b, 3);
    T[2][2] = sgn(c) * delta(c, 3);

    IJK begin2 = donorBeginIn;
    IJK delta1 = {recvEnd[0] - recvBegin[0], recvEnd[1] - recvBegin[1], recvEnd[2] - recvBegin[2]};

    auto apply = [&](const IJK& start) {
        IJK out{};
        for (int row = 0; row < 3; ++row) {
            long long sum = 0;
            for (int col = 0; col < 3; ++col) {
                sum += static_cast<long long>(T[row][col]) * delta1[col];
            }
            out[row] = sum + start[row];
        }
        return out;
    };

    IJK index2 = apply(begin2);
    bool redo = false;

    for (int dim = 0; dim < 3; ++dim) {
        if (index2[dim] < begin2[dim]) {
            begin2[dim] += 1;
            redo = true;
        }
    }

    if (redo) {
        index2 = apply(begin2);
    }

    donorBegin = begin2;
    donorEnd = index2;
}

} // namespace fs::conn
