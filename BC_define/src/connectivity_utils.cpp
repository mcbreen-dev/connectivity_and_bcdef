//─────────────────────────────────────────────────────────────
// File: src/connectivity_utils.cpp
// Shared math and indexing utilities for connectivity detection
//─────────────────────────────────────────────────────────────
//
// This file collects small, pure helper routines used by the structured
// multizone connectivity pipeline.
//
// Call sites (by file):
//   - connectivity_face_build.cpp
//       * plane_index(), volume_index(), face_dims(), face_cell_indices()
//       * face_grid_index() for FaceGrid population
//   - connectivity_hash.cpp
//       * expand_bounds(), bounds_overlap(), CellKey hashing helpers
//   - connectivity_verify.cpp
//       * face_dir_from_id(), is_i_face()/is_j_face()/is_k_face()
//       * signed_axis(), set_normal_transform(), fill_point_range()
//       * compute_donor_range()
//       * distance_sq() for tolerance tests
//   - connectivity_core.cpp
//       * is_left_handed(), is_left_handed_plot3d() validation gate
//
// The codebase uses 1-based structured indices (I,J,K) when referring to
// mesh vertices in CGNS/Plot3D conventions. Internal flattened arrays are
// 0-based, and the indexing helpers below perform the conversion.
#include "connectivity_internal.hpp"

#include <algorithm>
#include <cmath>

extern "C" {
    #include <cgnslib.h>
}

namespace fs::conn {

/*=====================================================================
// plane_index
//
// Convert 1-based structured coordinates (a,b) on a 2-D plane into a
// 0-based linear index into FacePlane::x/y/z.
//
// Inputs:
//   plane.dim1 : number of points in the fast-varying "a" direction
//   a, b       : 1-based plane coordinates (a ∈ [1..dim1], b ∈ [1..dim2])
//
// Output:
//   0-based index into FacePlane coordinate arrays.
//
// Layout:
//   index = (b-1)*dim1 + (a-1)
//=====================================================================*/
size_t plane_index(const FacePlane& plane, long long a, long long b)
{
    return static_cast<size_t>((b - 1) * plane.dim1 + (a - 1));
}

/*=====================================================================
// face_grid_index
//
// Convert 0-based face-local cell coordinates (u,v) into a 0-based linear
// index into FaceGrid::index and FaceGrid::consumed.
//
// Inputs:
//   grid.uLen : number of cells in u direction
//   u, v      : 0-based indices (u ∈ [0..uLen-1], v ∈ [0..vLen-1])
//
// Output:
//   0-based index into grid arrays.
//
// Layout:
//   idx = v*uLen + u
//=====================================================================*/
size_t face_grid_index(const FaceGrid& grid, int u, int v)
{
    return static_cast<size_t>(v) * grid.uLen + u;
}

/*=====================================================================
// volume_index
//
// Convert 1-based structured vertex indices (i,j,k) into a 0-based
// linear index into VolumeCoords::x/y/z.
//
// Inputs:
//   vol.ni, vol.nj : vertex dimensions (Ni,Nj) used for stride
//   i,j,k          : 1-based vertex indices
//
// Output:
//   0-based index into volume coordinate arrays.
//
// Layout (I fastest, then J, then K):
//   idx = (k-1)*Ni*Nj + (j-1)*Ni + (i-1)
//=====================================================================*/
size_t volume_index(const VolumeCoords& vol, long long i, long long j, long long k)
{
    return static_cast<size_t>((k - 1) * vol.ni * vol.nj + (j - 1) * vol.ni + (i - 1));
}

/*=====================================================================
// plot3d_index
//
// Same linearization as volume_index(), but using Plot3DZone accessors.
//
// Plot3DZone stores its vertex coordinates in three flat arrays x/y/z,
// with I fastest, then J, then K.
//
// Inputs:
//   zone.ni(), zone.nj() : used for stride
//   i,j,k                : 1-based vertex indices
//
// Output:
//   0-based index into zone.x/y/z.
//=====================================================================*/
size_t plot3d_index(const Plot3DZone& zone, long long i, long long j, long long k)
{
    return static_cast<size_t>((k - 1) * zone.ni() * zone.nj() + (j - 1) * zone.ni() + (i - 1));
}

/*=====================================================================
// is_left_handed
//
// Detect whether a structured CGNS zone has left-handed orientation.
//
// This routine is called from connectivity_core.cpp before connectivity
// generation for 3-D zones. A negative determinant implies that the
// local basis (di,dj,dk) at (1,1,1) is left-handed.
//
// Inputs:
//   fn, B : CGNS file and base IDs
//   z     : zone metadata (z.idx used for cg_coord_read)
//
// Data access:
//   Reads single-point RealDouble values from the CoordinateX/Y/Z arrays.
//   Each cg_coord_read reads exactly one double because start==end.
//
// 2D behavior:
//   If z.nk() == 1, dk is evaluated at k=1 so the determinant uses a
//   degenerate dk. connectivity_core.cpp skips this check for nk==1 zones.
//=====================================================================*/
bool is_left_handed(int fn, int B, const Zone& z)
{
    static constexpr const char* cname[3] =
        {"CoordinateX", "CoordinateY", "CoordinateZ"};

    // Read one vertex coordinate triple (X,Y,Z) at a specific (i,j,k).
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

    // Reference vertex and basis step vectors.
    Point p000 = node(1, 1, 1);
    Point di   = node(2, 1, 1);
    Point dj   = node(1, 2, 1);
    Point dk   = node(1, 1, (z.nk() > 1 ? 2 : 1));

    di.x -= p000.x; di.y -= p000.y; di.z -= p000.z;
    dj.x -= p000.x; dj.y -= p000.y; dj.z -= p000.z;
    dk.x -= p000.x; dk.y -= p000.y; dk.z -= p000.z;

    // det([di dj dk]) < 0 => left-handed orientation.
    double det = di.x * (dj.y * dk.z - dj.z * dk.y)
               - di.y * (dj.x * dk.z - dj.z * dk.x)
               + di.z * (dj.x * dk.y - dj.y * dk.x);
    return det < 0.0;
}

/*=====================================================================
// is_left_handed_plot3d
//
// Plot3D equivalent of is_left_handed(), using Plot3DZone::x/y/z arrays.
//
// Inputs:
//   z : Plot3D zone (1-based indexing convention for node() arguments)
//
// Output:
//   true if det([di dj dk]) < 0.
//
// Call site:
//   connectivity_core.cpp checks 3-D Plot3D zones and throws if left-handed.
//=====================================================================*/
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

/*=====================================================================
// distance_sq
//
// Squared Euclidean distance, used for tolerance comparisons without
// an expensive sqrt.
//
// Call sites:
//   - connectivity_verify.cpp (find_match_in_face and final checks)
//
// Inputs:
//   a,b : Point values with double x/y/z
// Output:
//   (ax-bx)^2 + (ay-by)^2 + (az-bz)^2
//=====================================================================*/
double distance_sq(const Point& a, const Point& b)
{
    double dx = a.x - b.x;
    double dy = a.y - b.y;
    double dz = a.z - b.z;
    return dx * dx + dy * dy + dz * dz;
}

/* Overload for raw coordinates to avoid temporary Point construction. */
double distance_sq(double ax, double ay, double az,
                   double bx, double by, double bz)
{
    double dx = ax - bx;
    double dy = ay - by;
    double dz = az - bz;
    return dx * dx + dy * dy + dz * dz;
}

/*=====================================================================
// expand_bounds
//
// Expand an axis-aligned bounding box (FaceBounds) to include a point.
// FaceBounds::valid is used to detect uninitialized bounds.
//
// Call sites:
//   - connectivity_hash.cpp builds per-face and global bounds
//=====================================================================*/
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

/*=====================================================================
// bounds_overlap
//
// Conservative AABB overlap test with tolerance expansion.
// Returns true if boxes overlap (or nearly overlap) in all axes.
//
// Call site:
//   - connectivity_hash.cpp final filter before generating FacePair list
//=====================================================================*/
bool bounds_overlap(const FaceBounds& a, const FaceBounds& b, double tol)
{
    if (!a.valid || !b.valid)
        return false;
    if (a.min.x > b.max.x + tol || b.min.x > a.max.x + tol) return false;
    if (a.min.y > b.max.y + tol || b.min.y > a.max.y + tol) return false;
    if (a.min.z > b.max.z + tol || b.min.z > a.max.z + tol) return false;
    return true;
}

/*=====================================================================
// cross / dot / norm / normalize
//
// Minimal vector math used in face summary construction and diagnostics.
//
// normalize() returns (0,0,0) for zero-length inputs.
//=====================================================================*/
Point cross(const Point& a, const Point& b)
{
    return Point{
        a.y * b.z - a.z * b.y,
        a.z * b.x - a.x * b.z,
        a.x * b.y - a.y * b.x
    };
}

double dot(const Point& a, const Point& b)
{
    return a.x * b.x + a.y * b.y + a.z * b.z;
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

/*=====================================================================
// CellKeyHash
//
// Hash functor for CellKey{ix,iy,iz} used as unordered_map key.
//
// Call sites:
//   - connectivity_hash.cpp (face_buckets and per-face FaceHash tables)
//   - connectivity_verify.cpp (FaceHash lookup)
//
// The combine scheme mixes ix, then incorporates iy and iz with a
// 64-bit constant and shift/xor steps.
//=====================================================================*/
size_t CellKeyHash::operator()(const CellKey& key) const
{
    size_t h = static_cast<size_t>(key.ix);
    h ^= static_cast<size_t>(key.iy + 0x9e3779b97f4a7c15ULL) + (h << 6) + (h >> 2);
    h ^= static_cast<size_t>(key.iz + 0x9e3779b97f4a7c15ULL) + (h << 6) + (h >> 2);
    return h;
}

/*=====================================================================
// quantize_coord / quantize_coord_offset
//
// Convert a coordinate into an integer cell index using floor().
//
// quantize_coord(v, cell):
//   returns floor(v / cell)
//
// quantize_coord_offset(v, origin, cell):
//   returns floor((v - origin) / cell)
//
// Call sites:
//   - connectivity_hash.cpp (face-level bucket index ranges)
//   - connectivity_hash.cpp / connectivity_verify.cpp via cell_key()
//=====================================================================*/
std::int64_t quantize_coord(double v, double cell)
{
    return static_cast<std::int64_t>(std::floor(v / cell));
}

std::int64_t quantize_coord_offset(double v, double origin, double cell)
{
    return static_cast<std::int64_t>(std::floor((v - origin) / cell));
}

/*=====================================================================
// cell_key
//
// Compute a 3-D integer cell key by quantizing coordinates.
// Two overloads are provided for Point and raw x/y/z.
//
// Call sites:
//   - connectivity_hash.cpp for per-face point hash construction
//   - connectivity_verify.cpp for match probe key generation
//=====================================================================*/
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

/*=====================================================================
// face_dir_from_id
//
// Map the connectivity pipeline's numeric face ID (0..5) into the public
// FaceDir enum used by BoundaryRecord / ConnPatch.
//
// Face numbering used throughout FaceBuildResult:
//   0=iMin, 1=iMax, 2=jMin, 3=jMax, 4=kMin, 5=kMax
//
// Call sites:
//   - connectivity_verify.cpp when emitting BoundaryRecord and ConnPatch
//=====================================================================*/
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

/* Face classification helpers used when building transforms. */
bool is_i_face(int face) { return face == 0 || face == 1; }
bool is_j_face(int face) { return face == 2 || face == 3; }
bool is_k_face(int face) { return face == 4 || face == 5; }

/*=====================================================================
// face_dims (Zone / Plot3DZone)
//
// Return the face-local cell-grid dimensions (uLen,vLen) for a given
// face of a structured zone.
//
// Definitions:
//   - uLen and vLen are counts of *cells* along the two varying axes.
//   - For a structured zone with vertex sizes (Ni,Nj,Nk):
//       cells are (Ni-1, Nj-1, Nk-1) in each axis.
//
// Face mapping (3D):
//   face 0/1 (iMin/iMax): u=J, v=K  => (Nj-1, Nk-1)
//   face 2/3 (jMin/jMax): u=I, v=K  => (Ni-1, Nk-1)
//   face 4/5 (kMin/kMax): u=I, v=J  => (Ni-1, Nj-1)
//
// 2D special-case (Nk==1):
//   Only iMin/iMax and jMin/jMax are considered; vLen is forced to 1 so
//   the rest of the pipeline can treat 2D faces as 1-row grids.
//   kMin/kMax return (0,0) and are skipped by callers.
//=====================================================================*/
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

/*=====================================================================
// face_cell_indices (Zone / Plot3DZone)
//
// Convert a face-local cell coordinate (u,v) to the corresponding
// 1-based vertex index (i,j,k) on that zone face.
//
// Inputs:
//   face : 0..5 in the same numbering as face_dir_from_id()
//   u,v  : 0-based cell coordinates on the face grid
//
// Outputs:
//   i,j,k : 1-based *vertex* indices, placed on the face plane.
//           The indices represent the "upper-right" vertex of the cell
//           in that direction, consistent with the build pipeline's
//           convention of storing face grid points at vertex locations.
//
// 2D behavior:
//   For Nk==1, k is pinned to 1 and v is ignored.
//=====================================================================*/
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

/*=====================================================================
// should_write
//
// Deterministic tie-break rule for emitting a single 1-to-1 connection
// for an interface that is encountered twice (once from each side).
//
// Call site:
//   - resolve_matches_*() in connectivity_verify.cpp
//
// Rule:
//   - Prefer smaller recvZone.
//   - If zones equal (should not happen for inter-zone links), prefer
//     smaller recvFace.
//=====================================================================*/
bool should_write(int recvZone, int donorZone, int recvFace, int donorFace)
{
    if (recvZone < donorZone) return true;
    if (recvZone > donorZone) return false;
    return recvFace < donorFace;
}

/*=====================================================================
// signed_axis
//
// Utility to express an axis as ±{1,2,3} depending on delta sign.
// axis is the unsigned axis id (1=I, 2=J, 3=K).
//
// Call site:
//   - resolve_matches_*() uses this while discovering donor index step
//     directions from adjacent points.
//=====================================================================*/
int signed_axis(int axis, int delta)
{
    return (delta >= 0) ? axis : -axis;
}

/*
CGNS-SIDS (GridConnectivity1to1_t::Transform):
Transform is shorthand for a signed-permutation (orthonormal) matrix T relating interface indices:
    IndexDonor = T * (IndexRecv - BeginRecv) + BeginDonor
Elements are signed integers with unique magnitudes in [1..IndexDimension]; in 3-D values may be 0, ±1, ±2, ±3.
For an interface plane, one component may be set to 0; it can be reconstructed from
the receiver/donor faces (min/max) and grid sizes if continuation off the interface is needed.

CGNS User Guide. CGNS/SIDS - Standard Interface Data Structures/8.2. 1-to-1 Interface Connectivity Structure Definition
https://cgns.org/standard/SIDS/multizone.html
*/

/*=====================================================================
// set_normal_transform
//
// Fill the transform component corresponding to the receiver's collapsed
// axis (the face normal direction), based on receiver face and donor face.
//
// This function writes exactly one of transform[0..2]; callers set the
// two tangential components before calling.
//
// Inputs:
//   transform : partially filled signed-axis mapping (size 3)
//   recvFace  : receiver face id (0..5)
//   donorFace : donor face id (0..5)
//
// Output:
//   transform[axis_of_recv_face] is set to ±donorAxis, where donorAxis is
//   1/2/3 for donor I/J/K, and the sign encodes whether the interface is
//   min-to-min, min-to-max, max-to-min, or max-to-max.
//
// Sign computation:
//   recvId = recvFace+1 and donorId = donorFace+1 gives IDs 1..6.
//   For the I/J/K min/max pairs, IDs differ by 1. The parity expression
//   produces -1 for same-side (min-min or max-max) and +1 for opposite
//   sides (min-max or max-min).
//
// 2D behavior:
//   In 2D paths, callers later force transform[2]=0 and force remaining
//   components to valid non-K axes.
//=====================================================================*/
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

/*=====================================================================
// fill_point_range (Zone / Plot3DZone)
//
// Convert a rectangular face-local cell patch [uStart..uEnd]×[vStart..vEnd]
// into CGNS-style vertex point-range begin/end (inclusive).
//
// Inputs:
//   face          : 0..5 face id
//   uStart/uEnd   : 0-based indices in u (cell grid coordinates)
//   vStart/vEnd   : 0-based indices in v (cell grid coordinates)
//
// Outputs:
//   begin/end     : IJK (long long[3]) with 1-based vertex indices.
//
// The "+1/+2" scheme:
//   - uStart is a cell index; the first vertex index on that edge is uStart+1.
//   - uEnd is the last cell; the far vertex index is (uEnd+1)+1 = uEnd+2.
//   This yields a vertex range that spans the cell patch.
//
// Face-specific mapping:
//   For each face, one axis is constant (the face plane), and the other
//   two axes vary using u0/u1 and v0/v1.
//
// 2D behavior:
//   For Nk==1, k is pinned to 1 and v maps to k is ignored by callers.
//=====================================================================*/
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

/*=====================================================================
// compute_donor_range
//
// Given:
//   - receiver begin/end (inclusive IJK range)
//   - a CGNS-style signed-axis transform (transform[3])
//   - a donor start point (donorBeginIn), which corresponds to recvBegin
//
// Compute the donor end point so that:
//   donorEnd = T * (recvEnd - recvBegin) + donorBegin
//
// transform encoding:
//   transform[d] ∈ {0, ±1, ±2, ±3}
//     magnitude chooses donor axis (1=I,2=J,3=K)
//     sign chooses direction along that axis
//   A value 0 indicates a collapsed dimension (2D case).
//
// Implementation steps:
//   1) Convert the transform vector into an explicit 3×3 matrix T
//      with entries in {-1,0,+1}.
//   2) Compute delta1 = recvEnd - recvBegin (inclusive range delta).
//   3) Apply T to delta1 and add donorBeginIn to obtain donorEnd.
//   4) If any donorEnd component is less than donorBegin, increment that
//      donorBegin component by 1 and recompute donorEnd.
//      This is a corrective step used by resolve_matches_*() because the
//      donor "startOther" point comes from a particular boundary vertex
//      and the patch begin may need to be shifted to keep begin<=end.
//
// Outputs:
//   donorBegin : adjusted donor begin (possibly shifted +1 on some axes)
//   donorEnd   : computed donor end corresponding to recvEnd
//
// Call sites:
//   - resolve_matches_cgns()
//   - resolve_matches_plot3d()
//=====================================================================*/
void compute_donor_range(const int transform[3], const IJK& recvBegin, const IJK& recvEnd,
                         const IJK& donorBeginIn, IJK& donorBegin, IJK& donorEnd)
{
    auto sgn = [](int x) { return x >= 0 ? 1 : -1; };
    auto delta = [](int x, int axis) { return (std::abs(x) == axis) ? 1 : 0; };

    int a = transform[0];
    int b = transform[1];
    int c = transform[2];

    // Expand the signed-permutation transform vector into a 3x3 matrix.
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

    // Delta of an inclusive receiver range in index space.
    IJK delta1 = {recvEnd[0] - recvBegin[0], recvEnd[1] - recvBegin[1], recvEnd[2] - recvBegin[2]};

    // Apply donorIndex = T*delta + start for the end corner.
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

    // Ensure donorBegin <= donorEnd component-wise by shifting begin and
    // recomputing the end if needed.
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
