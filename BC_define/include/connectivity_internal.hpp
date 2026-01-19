/*
  File: include/connectivity_internal.hpp

  Internal implementation API for connectivity detection.

  This header is included by:
    - src/connectivity_core.cpp
    - src/connectivity_face_build.cpp
    - src/connectivity_hash.cpp
    - src/connectivity_utils.cpp
    - src/connectivity_verify.cpp

  It defines the data structures and functions used to:
    1) enumerate structured boundary faces (CGNS + Plot3D)
    2) compute per-face-cell geometric "center" points (Point)
    3) build spatial hashes for fast proximity lookup of matching face cells
    4) collect point-to-point matches between candidate face pairs
    5) resolve those matches into:
         - ConnPatch entries (for GridConnectivity1to1 writing)
         - BoundaryRecord entries (for later BC assignment)

  Indexing and size conventions used throughout:
    - Zones/blocks:
        * build_faces_* uses 0-based "block" indices internally:
            block = [0 .. nBlocks-1]
        * external user-facing zone ids are 1-based (CGNS Zone::idx, Plot3DZone::idx).
        * ConnPatch stores recvZone/donorZone as 1-based.

    - Faces:
        * each block contributes up to 6 faces (IMIN, IMAX, JMIN, JMAX, KMIN, KMAX)
        * internal "face_id" = block*6 + face where face is 0..5
        * face_dir_from_id(face) converts 0..5 -> FaceDir

    - 2D handling:
        * 2D is represented by nk()==1
        * only IMIN/IMAX/JMIN/JMAX produce cells; K faces are treated as empty
        * some transforms collapse K (tvec[2]=0) when either side is 2D

    - (u,v) on a face:
        * each face has a local 2D grid with dimensions (uLen, vLen) measured in
          *cell* counts along the face:
            - uLen = (varying-axis length - 1)
            - vLen = (varying-axis length - 1) in 3D, or 1 in 2D
        * build stores u and v as int for each face cell center.

    - coordinate storage:
        * coordinates are stored as double in FacePlane / VolumeCoords / FaceData
        * memory usage scales with:
            - face planes: O(face vertex count) per face read
            - volume coords: O(ni*nj*nk) per zone if enabled by k_values thresholds
            - FaceData: O(total face cell count) across all faces

  CGNS coordinate reads:
    - FacePlane reads CoordinateX/Y/Z for a single boundary plane using cg_coord_read.
    - VolumeCoords reads CoordinateX/Y/Z for the full zone volume using cg_coord_read.
    - build_faces_cgns chooses between the two per zone using k-values and MEM_CAP_BYTES.
*/
#pragma once

#include "common.hpp"
#include "connectivity.hpp"
#include "mesh_io.hpp"
#include "plot3d_io.hpp"

#include <cstddef>
#include <cstdint>
#include <unordered_map>
#include <vector>

namespace bcdef::conn {


/*==============================================================================
  Small geometry primitives
==============================================================================*/

/*
  Point

  3D coordinate in double precision.

  Used for:
    - storing face-cell center coordinates (FaceData::x/y/z)
    - accumulating face centroid sums (FaceBuildResult::face_sum)
    - bounding boxes (FaceBounds::min/max)
    - normals (FaceSummary::normal)

  Size:
    - 3 * 8 bytes = 24 bytes on typical platforms
*/
struct Point {
    double x = 0.0;
    double y = 0.0;
    double z = 0.0;
};

/*
  FaceBounds

  Axis-aligned bounding box (AABB) for a face's cell-center points.

  Used for:
    - coarse overlap filtering between faces (bounds_overlap)
    - building face "bucket" candidates in build_hash_candidates()

  valid:
    - false until expand_bounds() is called the first time
*/
struct FaceBounds {
    bool valid = false;
    Point min{};
    Point max{};
};

/*==============================================================================
  Face-local coordinate containers (CGNS reading helpers)
==============================================================================*/

/*
  FacePlane

  Stores CoordinateX/Y/Z for one CGNS face plane (constant i, j, or k).

  Filled by:
    - read_face_plane(fn, B, zone, face) in src/connectivity_face_build.cpp

  face:
    - 0..5 (IMIN..KMAX)

  dim1, dim2:
    - vertex dimensions on that plane
      Example: IMIN plane has dim1 = Nj, dim2 = Nk

  x/y/z:
    - dim1*dim2 vertex values, laid out in plane_index() order:
        index = (b-1)*dim1 + (a-1)
      where (a,b) are the two varying indices on the plane.

  Sizes:
    - x/y/z each store (dim1*dim2) doubles.
*/
struct FacePlane {
    int face = -1;
    long long dim1 = 0;
    long long dim2 = 0;
    std::vector<double> x, y, z;
};

/*
  VolumeCoords

  Stores CoordinateX/Y/Z for the full structured zone.

  Filled by:
    - read_volume_coords(fn, B, zone) in src/connectivity_face_build.cpp

  ni/nj/nk:
    - vertex dimensions (Zone::ni/nj/nk)

  x/y/z:
    - ni*nj*nk doubles, indexed by volume_index():
        (k-1)*ni*nj + (j-1)*ni + (i-1)

  Used when build_faces_cgns decides "use_volume" for a zone.
*/
struct VolumeCoords {
    long long ni = 0;
    long long nj = 0;
    long long nk = 0;
    std::vector<double> x, y, z;
};

/*==============================================================================
  Packed face-cell point data (center points for matching)
==============================================================================*/

/*
  FaceData

  Structure-of-arrays (SoA) container for all face-cell centers across the mesh.

  Produced by:
    - FaceData::push(...) calls in build_faces_cgns() / build_faces_plot3d()

  Used by:
    - hashing (build_hash_candidates / find_match_in_face)
    - match verification and patch resolution (resolve_matches_*)

  Layout:
    - Each face cell center is one "record" at index idx.
    - The vectors store per-record attributes:
        block:  0-based block index (int)
        face:   0..5 face id within block (int)
        u,v:    face-local cell indices (int)
        i,j,k:  zone vertex index used as a reference anchor for donor-range derivation (int)
        x,y,z:  world coordinates of the center (double)

  Size per record (typical):
    - 7 ints (28 bytes if 4-byte int) + 3 doubles (24 bytes) = 52 bytes
      (plus vector overhead / alignment, not counted here)
*/
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

    /*
      reserve(n)
      Pre-allocates capacity in each vector to avoid reallocation.

      Called by:
        - build_faces_cgns() and build_faces_plot3d() before pushing centers.
    */
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

    /*
      size()
      Returns number of records stored (all vectors are kept the same length).
    */
    size_t size() const { return x.size(); }
    bool empty() const { return x.empty(); }

    /*
      push(...)
      Appends one face-cell center record.

      Returns:
        - index of inserted record (int), used as the canonical "cell id"
          throughout matching and resolution.
    */
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

/*==============================================================================
  Face grid / build products
==============================================================================*/

/*
  FaceGrid

  Dense index map from (u,v) -> FaceData record index for one face.

  Built in:
    - build_faces_cgns()
    - build_faces_plot3d()

  uLen, vLen:
    - dimensions of this face's cell grid

  index:
    - length uLen*vLen
    - index[face_grid_index(grid,u,v)] gives record id or -1

  consumed:
    - length uLen*vLen
    - used by resolve_matches_* to flood/scan rectangles and avoid re-processing.
    - values are 0/1 (bytes) for faster writes than bool bit-packing.
*/
struct FaceGrid {
    int uLen = 0;
    int vLen = 0;
    std::vector<int> index;
    std::vector<unsigned char> consumed;
};

/*
  FaceBuildResult

  Output of face enumeration and center computation.

  Produced by:
    - build_faces_cgns(mesh, show_progress)
    - build_faces_plot3d(zones)

  grids:
    - per-face dense index map; size = nBlocks*6

  bounds:
    - per-face AABB of cell centers; size = nBlocks*6

  cell_start / cell_count:
    - per-face slice into face_cells:
        * start = cell_start[face_id]
        * count = cell_count[face_id]
        * face_cells[start..start+count) are FaceData indices on that face

  face_cells:
    - concatenation of FaceData indices grouped per face

  faces:
    - FaceData storage of all centers

  face_sum / face_counts:
    - per-face accumulation of center sums and counts for centroid computation

  faces_loaded:
    - counter used for progress logging:
        * incremented by 6 for 3D zones, 4 for 2D zones
*/
struct FaceBuildResult {
    std::vector<FaceGrid> grids;
    std::vector<FaceBounds> bounds;

    std::vector<size_t> cell_start;     // size = face_count
    std::vector<int> cell_count;        // size = face_count
    std::vector<int> face_cells;        // length = total face cell count

    FaceData faces;

    std::vector<Point> face_sum;        // size = face_count
    std::vector<int> face_counts;       // size = face_count

    size_t faces_loaded = 0;
};

/*
  FaceSummary

  Per-face metadata derived from FaceBuildResult.

  Produced by:
    - build_face_summaries(build)

  Used by:
    - build_hash_candidates() (bounds for overlap and bucket building)
    - match collection (uLen/vLen flags for uv mapping)
    - logging/diagnostics (centroid/normal exists for future checks)

  normal_valid:
    - true only if uLen>1 and vLen>1 and a normal length > 0
      computed from three corner centers on the face.

  Sizes:
    - centroid/normal are Point (3 doubles each)
    - bounds is FaceBounds (bool + 2 Points)
*/
struct FaceSummary {
    int face_id = 0;        // block*6 + face
    int block = 0;          // 0-based
    int face = 0;           // 0..5
    int uLen = 0;
    int vLen = 0;
    bool normal_valid = false;
    Point normal{};
    Point centroid{};
    FaceBounds bounds{};
};

/*==============================================================================
  Spatial hashing for point matching
==============================================================================*/

/*
  CellKey

  Integer coordinates of a hash cell in a 3D uniform grid.

  Used as unordered_map key for:
    - FaceHash::cells (per-face: maps cell -> list of FaceData indices)
    - build_hash_candidates face_buckets (maps bucket -> list of face_ids)

  ix/iy/iz:
    - std::int64_t to support large coordinate magnitudes and bucket indices
*/
struct CellKey {
    std::int64_t ix = 0;
    std::int64_t iy = 0;
    std::int64_t iz = 0;

    bool operator==(const CellKey& other) const
    {
        return ix == other.ix && iy == other.iy && iz == other.iz;
    }
};

/*
  CellKeyHash
  Hash functor for CellKey.
  Implemented in src/connectivity_utils.cpp.
*/
struct CellKeyHash {
    size_t operator()(const CellKey& key) const;
};


/*
  FaceHash

  Hash table for one face's cell-center points.

  cells:
    - unordered_map<CellKey, vector<int>>
    - value list stores FaceData indices belonging to that cell

  Built by:
    - build_hash_candidates()
*/
struct FaceHash {
    std::unordered_map<CellKey, std::vector<int>, CellKeyHash> cells;
};

/*
  FacePair

  Candidate face pair to consider for matching.

  a,b:
    - face_id values (block*6 + face)
    - in build_hash_candidates, pairs are normalized to a<b before dedup.
*/
struct FacePair {
    int a = 0;
    int b = 0;
};


/*
  HashBuildResult

  Output of build_hash_candidates().

  hash_cell:
    - spatial cell size used for FaceHash construction
    - set to tol if tol>0, else 1.0

  face_hashes:
    - vector length face_count, per face hash of its points

  pairs:
    - candidate face pairs after global bucket filtering and bounds overlap test
*/
struct HashBuildResult {
    std::vector<FaceHash> face_hashes;
    std::vector<FacePair> pairs;
    double hash_cell = 1.0;
};

/*==============================================================================
  Matching and resolution helper records
==============================================================================*/

/*
  UVMatch

  Small sample record used by derive_uv_mapping().

  au,av:
    - (u,v) indices on face A

  bu,bv:
    - (u,v) indices on face B

  a_idx, b_idx:
    - FaceData record indices for the matched points

  Collected in:
    - collect_match_pairs() during "fast mapping" attempt
*/
struct UVMatch {
    int au = 0;
    int av = 0;
    int bu = 0;
    int bv = 0;
    int a_idx = -1;
    int b_idx = -1;
};

/*
  MatchPair

  Final point-to-point match record between two faces.

  a,b:
    - FaceData indices (record ids)
    - a is from face A and b is from face B for the processed candidate pair
    - collect_match_pairs() sorts output by (a,b)
*/
struct MatchPair {
    int a = 0;
    int b = 0;
};

/*
  NeighborInfo

  Stores the "other side" mapping for one FaceData record, after matches
  are filtered to be one-to-one.

  block:
    - 0-based donor block index
    - -1 means "no neighbor" (boundary)

  face:
    - 0..5 donor face id within donor block

  i,j,k:
    - donor zone indices recorded from FaceData::i/j/k
    - used as donorStart anchor in compute_donor_range()

  Built in:
    - resolve_matches_cgns()
    - resolve_matches_plot3d()

  Size:
    - 5 ints (typically 20 bytes)
*/
struct NeighborInfo {
    int block = -1;
    int face = -1;
    int i = -1;
    int j = -1;
    int k = -1;
};

/*
  BoundaryRecord

  Internal boundary output used by connectivity_core.cpp to form BoundaryPatch.

  zone:
    - 1-based zone index

  face:
    - FaceDir for the boundary face

  begin/end:
    - vertex-space IJK ranges (1-based), inclusive endpoints (as stored in IJK)
    - for partial patches: these enclose the patch rectangle in vertex indices
    - for full faces: these span the full face in vertex indices

  full:
    - true if the rectangle covers u:[0..uLen-1], v:[0..vLen-1]
      (the resolve code computes this exact condition)
*/
struct BoundaryRecord {
    int zone = 0;
    FaceDir face = FaceDir::IMIN;
    bool full = false;
    IJK begin{};
    IJK end{};
};

/*==============================================================================
  Core algorithms (implemented in .cpp files)
==============================================================================*/

/*------------------------ Indexing helpers ------------------------*/

/*
  plane_index(plane, a, b)

  Maps 1-based (a,b) indices on a FacePlane into a 0-based linear index.

  Used by:
    - face_center_from_plane() in src/connectivity_face_build.cpp
*/
size_t plane_index(const FacePlane& plane, long long a, long long b);

/*
  face_grid_index(grid, u, v)

  Maps 0-based (u,v) in a FaceGrid to a 0-based linear index.

  Used across:
    - build_faces_* (writing grid.index)
    - collect_match_pairs() (gridB lookup)
    - resolve_matches_* (scan rectangles; consumed marking)
*/
size_t face_grid_index(const FaceGrid& grid, int u, int v);

/*
  volume_index(vol, i, j, k)

  Maps 1-based (i,j,k) to 0-based linear index into VolumeCoords arrays.

  Used by:
    - face_center_from_volume()
*/
size_t volume_index(const VolumeCoords& vol, long long i, long long j, long long k);

/*
  plot3d_index(zone, i, j, k)

  Maps 1-based (i,j,k) to 0-based linear index into Plot3DZone x/y/z arrays.

  Used by:
    - face_center_from_plot3d()
    - is_left_handed_plot3d()
*/
size_t plot3d_index(const Plot3DZone& zone, long long i, long long j, long long k);


/*------------------------ Orientation / handedness ------------------------*/

/*
  is_left_handed(fn, B, zone)

  Reads 4 CGNS points:
    p000 = (1,1,1)
    di   = (2,1,1) - p000
    dj   = (1,2,1) - p000
    dk   = (1,1,2) - p000   (or (1,1,1) if nk==1)

  Computes determinant det = dot(di, cross(dj, dk)).
  Returns true if det < 0.

  Called by:
    - ConnectivityDetector::run_collect() (src/connectivity_core.cpp)
      for every 3D zone (nk>1). Left-handed zones throw.
*/
bool is_left_handed(int fn, int B, const Zone& z);

/*
  is_left_handed_plot3d(zone)

  Same test as is_left_handed, but using Plot3DZone arrays.

  Called by:
    - ConnectivityDetector::run_plot3d() (src/connectivity_core.cpp)
*/
bool is_left_handed_plot3d(const Plot3DZone& z);

/*------------------------ Small math utilities ------------------------*/
double distance_sq(const Point& a, const Point& b);
double distance_sq(double ax, double ay, double az, double bx, double by, double bz);
void expand_bounds(FaceBounds& bounds, const Point& p);
bool bounds_overlap(const FaceBounds& a, const FaceBounds& b, double tol);
Point cross(const Point& a, const Point& b);
double norm(const Point& a);
Point normalize(const Point& a);
double dot(const Point& a, const Point& b);

/*------------------------ Hash utilities ------------------------*/

/*
  quantize_coord(v, cell)

  Returns floor(v / cell) as int64.

  Used by:
    - cell_key(...)
*/
std::int64_t quantize_coord(double v, double cell);

/*
  quantize_coord_offset(v, origin, cell)

  Returns floor((v - origin) / cell) as int64.

  Used by:
    - build_hash_candidates() when bucketing faces by their bounds relative to global min.
*/
std::int64_t quantize_coord_offset(double v, double origin, double cell);

/*
  cell_key(p, cell) / cell_key(x,y,z,cell)

  Converts coordinates into a CellKey for a given cell size.
*/
CellKey cell_key(const Point& p, double cell);
CellKey cell_key(double x, double y, double z, double cell);

/*------------------------ Face topology helpers ------------------------*/

FaceDir face_dir_from_id(int face);

bool is_i_face(int face);
bool is_j_face(int face);
bool is_k_face(int face);

/*
  face_dims(z, face, uLen, vLen)

  Computes face-local cell dimensions for:
    - CGNS Zone overload
    - Plot3DZone overload

  Exact behavior matches src/connectivity_utils.cpp:
    2D (nk==1):
      IMIN/IMAX -> uLen = nj-1, vLen = 1
      JMIN/JMAX -> uLen = ni-1, vLen = 1
      KMIN/KMAX -> uLen = 0, vLen = 0
    3D:
      IMIN/IMAX -> uLen = nj-1, vLen = nk-1
      JMIN/JMAX -> uLen = ni-1, vLen = nk-1
      KMIN/KMAX -> uLen = ni-1, vLen = nj-1
*/
void face_dims(const Zone& z, int face, int& uLen, int& vLen);
void face_dims(const Plot3DZone& z, int face, int& uLen, int& vLen);

/*
  face_cell_indices(z, face, u, v, i, j, k)

  Computes a representative vertex index (i,j,k) for the face cell at (u,v).

  Used in:
    - build_faces_cgns()
    - build_faces_plot3d()

  These i/j/k values are not used for geometry. They are used later as:
    - NeighborInfo::i/j/k donorStart anchor in resolve_matches_*()

  The mapping is constant-face:
    face 0 (IMIN): i=1,        j=u+1, k=v+1 (or 1 in 2D)
    face 1 (IMAX): i=ni,       j=u+1, k=v+1
    face 2 (JMIN): j=1,        i=u+1, k=v+1
    face 3 (JMAX): j=nj,       i=u+1, k=v+1
    face 4 (KMIN): k=1,        i=u+1, j=v+1
    face 5 (KMAX): k=nk,       i=u+1, j=v+1
*/
void face_cell_indices(const Zone& z, int face, int u, int v, int& i, int& j, int& k);
void face_cell_indices(const Plot3DZone& z, int face, int u, int v, int& i, int& j, int& k);

/*------------------------ Connection writing policy ------------------------*/

/*
  should_write(recvZone, donorZone, recvFace, donorFace)

  Enforces a deterministic direction so each 1-to-1 connection is emitted once.

  Rule (src/connectivity_utils.cpp):
    - write if recvZone < donorZone
    - if zones equal, write if recvFace < donorFace

  Used in:
    - resolve_matches_cgns()
    - resolve_matches_plot3d()
*/
bool should_write(int recvZone, int donorZone, int recvFace, int donorFace);

/*
  signed_axis(axis, delta)

  Returns +axis if delta>=0 else -axis.
  axis is 1,2,3 for I,J,K.

  Used in:
    - resolve_matches_*() while detecting which donor coordinate changes along u/v.
*/
int signed_axis(int axis, int delta);

/*
  set_normal_transform(transform, recvFace, donorFace)

  Writes the normal component of the CGNS transform vector based on:
    - which axis is normal for donorFace
    - sign convention derived from (recvId + donorId) parity
      (recvId = recvFace+1, donorId = donorFace+1)

  Used in:
    - resolve_matches_*() after setting the two tangential mapping entries.
*/
void set_normal_transform(int transform[3], int recvFace, int donorFace);

/*------------------------ Range construction ------------------------*/

/*
  fill_point_range(z, face, uStart, vStart, uEnd, vEnd, begin, end)

  Converts a face-local cell rectangle [uStart..uEnd], [vStart..vEnd]
  into a vertex-space IJK PointRange.

  Key detail:
    - uStart/vStart are cell indices; the vertex range spans +1 in each direction.
    - implementation uses:
        u0 = uStart + 1
        u1 = uEnd   + 2
        v0 = vStart + 1
        v1 = vEnd   + 2

  begin/end are inclusive endpoints in vertex index space.

  Used in:
    - resolve_matches_*() to populate:
        * BoundaryRecord begin/end for boundary patches
        * ConnPatch recvRange for interior patches
*/
void fill_point_range(const Zone& z, int face, int uStart, int vStart, int uEnd, int vEnd,
    IJK& begin, IJK& end);
void fill_point_range(const Plot3DZone& z, int face, int uStart, int vStart, int uEnd, int vEnd,
    IJK& begin, IJK& end);

/*
  compute_donor_range(transform, recvBegin, recvEnd, donorBeginIn, donorBegin, donorEnd)

  Computes donor range endpoints from:
    - recvBegin/recvEnd (vertex-space)
    - donorBeginIn anchor (vertex-space)
    - transform vector (CGNS transform semantics, ±1..±3 and possibly 0 for 2D)

  Implementation details:
    - constructs 3x3 signed permutation matrix T from transform[3]
    - computes delta1 = recvEnd - recvBegin (per axis)
    - donorEnd = donorBeginIn + T * delta1 (with an adjustment step if any axis flips)

  Used in:
    - resolve_matches_*() to populate ConnPatch donorRange
*/
void compute_donor_range(const int transform[3],
                         const IJK& recvBegin, const IJK& recvEnd,
                         const IJK& donorBeginIn,
                         IJK& donorBegin, IJK& donorEnd);


/*==============================================================================
  Face extraction and summary building
==============================================================================*/

/*
  face_center_from_plane(plane, zone, face, u, v)

  Computes the geometric center of one face cell using FacePlane vertex samples.

  u and v are 0-based cell indices.
  For 3D faces, averages 4 vertices (quad).
  For 2D faces, averages 2 vertices (edge).
*/
Point face_center_from_plane(const FacePlane& plane, const Zone& z, int face, int u, int v);

/*
  read_face_plane(fn, B, zone, face)

  Reads CoordinateX/Y/Z for a single boundary face plane.

  Returns:
    - FacePlane with dim1/dim2 and x/y/z filled
    - if face has non-positive dimensions, returns an empty plane with face set.
*/
FacePlane read_face_plane(int fn, int B, const Zone& z, int face);

/*
  read_volume_coords(fn, B, zone)

  Reads CoordinateX/Y/Z for the full zone volume.
  Used only when build_faces_cgns selects volume-based center computation.
*/
VolumeCoords read_volume_coords(int fn, int B, const Zone& z);

/*
  face_center_from_volume(vol, zone, face, u, v)

  Computes the center of one face cell by sampling from VolumeCoords.

  This avoids reading 6 separate FacePlane arrays when volume reads are cheaper
  under the chosen heuristic (k-values and memory cap).
*/
Point face_center_from_volume(const VolumeCoords& vol, const Zone& z, int face, int u, int v);

/*
  face_center_from_plot3d(zone, face, u, v)

  Computes the center of one face cell from Plot3DZone x/y/z arrays.

  Used by:
    - build_faces_plot3d()
*/
Point face_center_from_plot3d(const Plot3DZone& zone, int face, int u, int v);

/*
  build_faces_cgns(mesh, show_progress)

  Builds FaceBuildResult for a CGNS mesh:
    - decides per-zone between reading full volume vs individual face planes
    - enumerates face cells (u,v) and stores their centers in FaceData
    - fills FaceGrid (u,v -> record id), bounds, and per-face slices

  Called by:
    - ConnectivityDetector::run_collect()
*/
FaceBuildResult build_faces_cgns(Mesh& mesh, bool show_progress);

/*
  build_faces_plot3d(zones)

  Builds FaceBuildResult for Plot3D zones using Plot3DZone arrays.

  Called by:
    - ConnectivityDetector::run_plot3d()
*/
FaceBuildResult build_faces_plot3d(const std::vector<Plot3DZone>& zones);

/*
  build_face_summaries(build)

  Computes:
    - FaceSummary bounds (copied from build.bounds)
    - centroid from face_sum/face_counts
    - an approximate normal from 3 corner centers when (uLen>1 && vLen>1)

  Called by:
    - ConnectivityDetector::run_collect()
    - ConnectivityDetector::run_plot3d()
*/
std::vector<FaceSummary> build_face_summaries(const FaceBuildResult& build);

/*==============================================================================
  Hash candidate generation and matching
==============================================================================*/

/*
  build_hash_candidates(build, summaries, tol)

  Produces:
    - per-face FaceHash tables keyed by cell_key(center, hash_cell)
    - candidate face pairs based on coarse bounds bucket overlap

  Two spatial grids are used:
    1) point hash grid per face:
        - cell size = hash_cell = (tol>0 ? tol : 1)
        - FaceHash stores lists of FaceData indices in each cell

    2) global face bucket grid:
        - bucket size = max(tol*4, global_extent/32), fallback to 1
        - each face inserts into all buckets its bounds overlaps
        - candidate pairs are generated per bucket, deduped, then filtered by bounds_overlap

  Called by:
    - ConnectivityDetector::run_collect()
    - ConnectivityDetector::run_plot3d()
*/
HashBuildResult build_hash_candidates(const FaceBuildResult& build,
                                      const std::vector<FaceSummary>& summaries,
                                      double tol);

/*
  collect_match_pairs(build, summaries, hash, tol, threads)

  For each candidate FacePair:
    - chooses smaller face as "A" for fewer probe points
    - attempts to derive an exact uv mapping using a few matched samples:
        derive_uv_mapping()
      If mapping is found:
        - maps every (u,v) in A to (bu,bv) in B and validates distance <= tol
      Otherwise:
        - performs per-point hash neighbor search (find_match_in_face)

  Returns:
    - vector<MatchPair> sorted by (a,b)

  Called by:
    - ConnectivityDetector::run_collect()
    - ConnectivityDetector::run_plot3d()
*/
std::vector<MatchPair> collect_match_pairs(const FaceBuildResult& build,
                                           const std::vector<FaceSummary>& summaries,
                                           const HashBuildResult& hash,
                                           double tol,
                                           int threads);

/*
  derive_uv_mapping(matches, need_u, need_v, ...)

  Attempts to infer a linear mapping between face A (u,v) and face B (u,v):
    [bu] = origin.bu + du*mu_u + dv*mv_u
    [bv] = origin.bv + du*mu_v + dv*mv_v

  The coefficients are constrained to axis-aligned ±1 swaps/flips:
    - mu_* and mv_* are in {-1,0,+1}
    - u and v must map to orthogonal directions

  Used in:
    - collect_match_pairs()
*/                         
bool derive_uv_mapping(const std::vector<UVMatch>& matches,
                       bool need_u,
                       bool need_v,
                       int& mu_u,
                       int& mu_v,
                       int& mv_u,
                       int& mv_v);

/*
  find_match_in_face(hash, faces, tx,ty,tz, tol2, cell)

  Finds the closest FaceData point in face B to the target point (tx,ty,tz):
    - probes the 27 neighboring hash cells around cell_key(tx,ty,tz,cell)
    - selects the candidate with minimum squared distance
    - returns candidate index if best distance <= tol2

  Used in:
    - collect_match_pairs() (fallback mode and seed sample gathering)
*/
int find_match_in_face(const FaceHash& hash,
                       const FaceData& faces,
                       double tx,
                       double ty,
                       double tz,
                       double tol2,
                       double cell);

/*==============================================================================
  Patch resolution (matches -> boundary records + ConnPatch)
==============================================================================*/

/*
  resolve_matches_cgns(mesh, log, build, matches, boundary_records, conn_patches)

  Converts point matches into:
    - boundary patches: BoundaryRecord (zone, face, begin/end, full)
    - interior patches: ConnPatch entries with explicit recv/donor ranges and transform

  Steps:
    1) Build otherInfo[] array (NeighborInfo per FaceData record).
       Each FaceData record is matched at most once:
         if either side already has a neighbor, the match is ignored.

    2) For each block and each face:
       - scan the (u,v) grid into maximal contiguous rectangles where NeighborInfo
         is either boundary (-1) or consistent donor (same block+face with unit steps)
       - mark consumed cells to avoid reprocessing
       - for each rectangle:
           * compute recvBegin/recvEnd (vertex-space) via fill_point_range()
           * if boundary: emit BoundaryRecord
           * if interior:
               - derive tangential transform components from observed donor steps
               - set normal component via set_normal_transform()
               - compute donorBegin/donorEnd via compute_donor_range()
               - apply should_write() to emit connection once
               - append ConnPatch

  Called by:
    - ConnectivityDetector::run_collect()
*/
void resolve_matches_cgns(const Mesh& mesh, Logger& log,
                          FaceBuildResult& build,
                          const std::vector<MatchPair>& matches,
                          std::vector<BoundaryRecord>& boundary_records,
                          std::vector<ConnPatch>& conn_patches);

/*
  resolve_matches_plot3d(zones, log, build, matches, boundary_records, conn_patches)

  Same algorithm as resolve_matches_cgns, but:
    - uses Plot3DZone dimension queries instead of CGNS Zone
    - donor zone references come from the Plot3D zone list

  Called by:
    - ConnectivityDetector::run_plot3d()
*/
void resolve_matches_plot3d(const std::vector<Plot3DZone>& zones, Logger& log,
                            FaceBuildResult& build,
                            const std::vector<MatchPair>& matches,
                            std::vector<BoundaryRecord>& boundary_records,
                            std::vector<ConnPatch>& conn_patches);

} // namespace bcdef::conn
