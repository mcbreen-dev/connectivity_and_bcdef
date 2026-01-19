//─────────────────────────────────────────────────────────────
// File: src/connectivity_face_build.cpp
//
// Face extraction and metadata construction.
//
// This file turns each *zone boundary face* into a dense 2D grid
// of “face cells” (topological quads in 3D, edges in 2D) and stores
// one geometric representative per face cell: its center point.
//
// The connectivity algorithm uses these centers for:
//   - per-face bounding boxes (FaceBounds) and centroids
//   - spatial hashing of points (connectivity_hash.cpp)
//   - point-to-point matching between candidate face pairs
//     (connectivity_verify.cpp)
//
// Output types produced here:
//   - FacePlane:     coordinates for a single constant-index plane
//   - VolumeCoords:  full zone coordinate volume (CoordinateX/Y/Z)
//   - FaceBuildResult:
//       * FaceData        build.faces: SoA arrays of all face-cell centers
//       * FaceGrid        build.grids: per-face (u,v)->index mapping into FaceData
//       * face_cells      flat list of FaceData indices for each face
//       * bounds/sums     per-face bounding boxes and centroid accumulators
//       * faces_loaded    progress counter (# of faces read from input)
//
// Coordinate index conventions:
//   - CGNS uses 1-based indices for i,j,k.
//   - Internal (u,v) indices in FaceGrid are 0-based.
//   - Helper functions convert (face,u,v) to CGNS (i,j,k) indices.
//
// 2D zones (nk == 1):
//   - Only faces {iMin, iMax, jMin, jMax} are meaningful.
//   - Faces kMin/kMax are treated as empty (uLen/vLen set to 0).
//   - vLen is set to 1 for 2D faces so the grid is still 2D in storage.
//─────────────────────────────────────────────────────────────
#include "connectivity_internal.hpp"
#include "k_values.hpp"

#include <algorithm>

extern "C" {
    #include <cgnslib.h>
}

namespace fs::conn {

/*=====================================================================
  face_center_from_plane

  Compute the geometric center of one face cell (u,v) from a FacePlane.

  Parameters:
    - plane: holds CoordinateX/Y/Z for a constant-index plane
             (i, j, or k fixed depending on face)
    - z:     zone metadata (ni/nj/nk); only used to detect 2D vs 3D
    - face:  integer face id:
              0=iMin, 1=iMax, 2=jMin, 3=jMax, 4=kMin, 5=kMax
    - u,v:   0-based cell coordinates in the face-local grid

  Returns:
    - Point (double x,y,z) which is the average of the vertices that bound
      the face cell:
        * 3D: average of 4 corner vertices of the quad
        * 2D: average of 2 endpoints of the edge (degenerate quad)
=====================================================================*/
Point face_center_from_plane(const FacePlane& plane, const Zone& z, int face, int u, int v)
{
    bool is2D = (z.nk() == 1);

    // Convert 0-based (u,v) to 1-based indices used by plane_index().
    // The plane stores vertex coordinates, so moving one cell in u/v
    // moves one vertex index.
    long long u0 = u + 1;
    long long v0 = v + 1;

    /*---------------------------------------------------------
      2D: face cells are edges, so each “cell center” is the
          midpoint of 2 vertices.
    ---------------------------------------------------------*/
    if (is2D) {
        if (face == 0 || face == 1) {
            // iMin / iMax: edge varies along j
            long long j = u0;
            size_t p1 = plane_index(plane, j, 1);
            size_t p2 = plane_index(plane, j + 1, 1);
            return Point{(plane.x[p1] + plane.x[p2]) * 0.5,
                         (plane.y[p1] + plane.y[p2]) * 0.5,
                         (plane.z[p1] + plane.z[p2]) * 0.5};
        }
        if (face == 2 || face == 3) {
            // jMin / jMax: edge varies along i
            long long i = u0;
            size_t p1 = plane_index(plane, i, 1);
            size_t p2 = plane_index(plane, i + 1, 1);
            return Point{(plane.x[p1] + plane.x[p2]) * 0.5,
                         (plane.y[p1] + plane.y[p2]) * 0.5,
                         (plane.z[p1] + plane.z[p2]) * 0.5};
        }
        // kMin/kMax are not used in 2D; return a placeholder.
        return Point{0, 0, 0};
    }

    /*---------------------------------------------------------
      3D: each face cell is a quad; average the 4 vertices.

      Face-local (u,v) axes depend on which face is fixed:
        - i faces: u=j, v=k
        - j faces: u=i, v=k
        - k faces: u=i, v=j
    ---------------------------------------------------------*/
    if (face == 0 || face == 1) {
        long long j = u0;
        long long k = v0;
        size_t p1 = plane_index(plane, j, k);
        size_t p2 = plane_index(plane, j + 1, k);
        size_t p3 = plane_index(plane, j, k + 1);
        size_t p4 = plane_index(plane, j + 1, k + 1);
        return Point{(plane.x[p1] + plane.x[p2] + plane.x[p3] + plane.x[p4]) * 0.25,
                     (plane.y[p1] + plane.y[p2] + plane.y[p3] + plane.y[p4]) * 0.25,
                     (plane.z[p1] + plane.z[p2] + plane.z[p3] + plane.z[p4]) * 0.25};
    }

    if (face == 2 || face == 3) {
        long long i = u0;
        long long k = v0;
        size_t p1 = plane_index(plane, i, k);
        size_t p2 = plane_index(plane, i + 1, k);
        size_t p3 = plane_index(plane, i, k + 1);
        size_t p4 = plane_index(plane, i + 1, k + 1);
        return Point{(plane.x[p1] + plane.x[p2] + plane.x[p3] + plane.x[p4]) * 0.25,
                     (plane.y[p1] + plane.y[p2] + plane.y[p3] + plane.y[p4]) * 0.25,
                     (plane.z[p1] + plane.z[p2] + plane.z[p3] + plane.z[p4]) * 0.25};
    }

    // k faces
    long long i = u0;
    long long j = v0;
    size_t p1 = plane_index(plane, i, j);
    size_t p2 = plane_index(plane, i + 1, j);
    size_t p3 = plane_index(plane, i, j + 1);
    size_t p4 = plane_index(plane, i + 1, j + 1);
    return Point{(plane.x[p1] + plane.x[p2] + plane.x[p3] + plane.x[p4]) * 0.25,
                 (plane.y[p1] + plane.y[p2] + plane.y[p3] + plane.y[p4]) * 0.25,
                 (plane.z[p1] + plane.z[p2] + plane.z[p3] + plane.z[p4]) * 0.25};
}

/*=====================================================================
  read_face_plane (CGNS)

  Reads CoordinateX/Y/Z for exactly one boundary plane of one zone.

  The plane is defined by "face" (0..5) and is represented as:
    - plane.dim1, plane.dim2: vertex dimensions of that plane
    - plane.x/y/z: vectors of length dim1*dim2
    - plane.face: stored for bookkeeping

  CGNS calls:
    cg_coord_read(..., s, e, plane.x.data()) for each coord component.

  Indexing:
    - s/e are 1-based inclusive start/end indices.
    - One of {i,j,k} is collapsed (s[axis]==e[axis]) to read a plane.
=====================================================================*/
FacePlane read_face_plane(int fn, int B, const Zone& z, int face)
{
    FacePlane plane;
    plane.face = face;

    // Full zone bounds, then collapse one axis depending on which face.
    cgsize_t s[3] = {1, 1, 1};
    cgsize_t e[3] = {static_cast<cgsize_t>(z.ni()),
                     static_cast<cgsize_t>(z.nj()),
                     static_cast<cgsize_t>(z.nk())};

    if (face == 0) { // iMin
        s[0] = 1;
        e[0] = 1;
        plane.dim1 = z.nj(); // j vertices
        plane.dim2 = z.nk(); // k vertices
    } else if (face == 1) { // iMax
        s[0] = static_cast<cgsize_t>(z.ni());
        e[0] = static_cast<cgsize_t>(z.ni());
        plane.dim1 = z.nj();
        plane.dim2 = z.nk();
    } else if (face == 2) { // jMin
        s[1] = 1;
        e[1] = 1;
        plane.dim1 = z.ni(); // i vertices
        plane.dim2 = z.nk(); // k vertices
    } else if (face == 3) { // jMax
        s[1] = static_cast<cgsize_t>(z.nj());
        e[1] = static_cast<cgsize_t>(z.nj());
        plane.dim1 = z.ni();
        plane.dim2 = z.nk();
    } else if (face == 4) { // kMin
        s[2] = 1;
        e[2] = 1;
        plane.dim1 = z.ni(); // i vertices
        plane.dim2 = z.nj(); // j vertices
    } else { // kMax
        s[2] = static_cast<cgsize_t>(z.nk());
        e[2] = static_cast<cgsize_t>(z.nk());
        plane.dim1 = z.ni();
        plane.dim2 = z.nj();
    }

    // A plane with non-positive dimensions is treated as empty.
    if (plane.dim1 <= 0 || plane.dim2 <= 0)
        return plane;

    // dim1*dim2 doubles per coordinate component.
    size_t count = static_cast<size_t>(plane.dim1 * plane.dim2);
    plane.x.resize(count);
    plane.y.resize(count);
    plane.z.resize(count);

    CG_CALL(cg_coord_read(fn, B, z.idx, "CoordinateX", CGNS_ENUMV(RealDouble),
                          s, e, plane.x.data()), "cg_coord_read X failed");
    CG_CALL(cg_coord_read(fn, B, z.idx, "CoordinateY", CGNS_ENUMV(RealDouble),
                          s, e, plane.y.data()), "cg_coord_read Y failed");
    CG_CALL(cg_coord_read(fn, B, z.idx, "CoordinateZ", CGNS_ENUMV(RealDouble),
                          s, e, plane.z.data()), "cg_coord_read Z failed");

    return plane;
}

/*=====================================================================
  read_volume_coords (CGNS)

  Reads the full CoordinateX/Y/Z volume for a zone into SoA arrays.

  VolumeCoords storage:
    - ni, nj, nk: long long zone dimensions (vertex counts)
    - x/y/z: vectors of length ni*nj*nk

  Memory sizing:
    - count = ni*nj*nk elements
    - bytes per coordinate = count * sizeof(double)
    - total bytes = 3 * count * sizeof(double)

  This path is selected by build_faces_cgns when the heuristic in
  k_values.hpp decides the full-volume read is cheaper than reading
  boundary planes face-by-face, and when the memory cap permits it.
=====================================================================*/
VolumeCoords read_volume_coords(int fn, int B, const Zone& z)
{
    VolumeCoords vol;
    vol.ni = z.ni();
    vol.nj = z.nj();
    vol.nk = z.nk();

    size_t count = static_cast<size_t>(vol.ni * vol.nj * vol.nk);
    vol.x.resize(count);
    vol.y.resize(count);
    vol.z.resize(count);

    cgsize_t s[3] = {1, 1, 1};
    cgsize_t e[3] = {static_cast<cgsize_t>(vol.ni),
                     static_cast<cgsize_t>(vol.nj),
                     static_cast<cgsize_t>(vol.nk)};

    CG_CALL(cg_coord_read(fn, B, z.idx, "CoordinateX", CGNS_ENUMV(RealDouble),
                          s, e, vol.x.data()), "cg_coord_read X failed");
    CG_CALL(cg_coord_read(fn, B, z.idx, "CoordinateY", CGNS_ENUMV(RealDouble),
                          s, e, vol.y.data()), "cg_coord_read Y failed");
    CG_CALL(cg_coord_read(fn, B, z.idx, "CoordinateZ", CGNS_ENUMV(RealDouble),
                          s, e, vol.z.data()), "cg_coord_read Z failed");

    return vol;
}

/*=====================================================================
  face_center_from_volume

  Same geometric output as face_center_from_plane, but the vertex
  coordinates come from VolumeCoords.

  VolumeCoords is indexed with volume_index(), which maps 1-based i,j,k
  to a 0-based flat index:
    idx = (k-1)*ni*nj + (j-1)*ni + (i-1)

  2D zones (nk==1) treat k as always 1.
=====================================================================*/
Point face_center_from_volume(const VolumeCoords& vol, const Zone& z, int face, int u, int v)
{
    bool is2D = (z.nk() == 1);
    long long u0 = u + 1;
    long long v0 = v + 1;

    auto node = [&](long long i, long long j, long long k) {
        size_t idx = volume_index(vol, i, j, k);
        return Point{vol.x[idx], vol.y[idx], vol.z[idx]};
    };

    if (is2D) {
        if (face == 0 || face == 1) {
            long long j = u0;
            Point p1 = node(face == 0 ? 1 : vol.ni, j, 1);
            Point p2 = node(face == 0 ? 1 : vol.ni, j + 1, 1);
            return Point{(p1.x + p2.x) * 0.5, (p1.y + p2.y) * 0.5, (p1.z + p2.z) * 0.5};
        }
        if (face == 2 || face == 3) {
            long long i = u0;
            Point p1 = node(i, face == 2 ? 1 : vol.nj, 1);
            Point p2 = node(i + 1, face == 2 ? 1 : vol.nj, 1);
            return Point{(p1.x + p2.x) * 0.5, (p1.y + p2.y) * 0.5, (p1.z + p2.z) * 0.5};
        }
        return Point{0, 0, 0};
    }

    if (face == 0 || face == 1) {
        long long j = u0;
        long long k = v0;
        long long i = (face == 0) ? 1 : vol.ni;
        Point p1 = node(i, j, k);
        Point p2 = node(i, j + 1, k);
        Point p3 = node(i, j, k + 1);
        Point p4 = node(i, j + 1, k + 1);
        return Point{(p1.x + p2.x + p3.x + p4.x) * 0.25,
                     (p1.y + p2.y + p3.y + p4.y) * 0.25,
                     (p1.z + p2.z + p3.z + p4.z) * 0.25};
    }

    if (face == 2 || face == 3) {
        long long i = u0;
        long long k = v0;
        long long j = (face == 2) ? 1 : vol.nj;
        Point p1 = node(i, j, k);
        Point p2 = node(i + 1, j, k);
        Point p3 = node(i, j, k + 1);
        Point p4 = node(i + 1, j, k + 1);
        return Point{(p1.x + p2.x + p3.x + p4.x) * 0.25,
                     (p1.y + p2.y + p3.y + p4.y) * 0.25,
                     (p1.z + p2.z + p3.z + p4.z) * 0.25};
    }

    long long i = u0;
    long long j = v0;
    long long k = (face == 4) ? 1 : vol.nk;
    Point p1 = node(i, j, k);
    Point p2 = node(i + 1, j, k);
    Point p3 = node(i, j + 1, k);
    Point p4 = node(i + 1, j + 1, k);
    return Point{(p1.x + p2.x + p3.x + p4.x) * 0.25,
                 (p1.y + p2.y + p3.y + p4.y) * 0.25,
                 (p1.z + p2.z + p3.z + p4.z) * 0.25};
}

/*=====================================================================
  face_center_from_plot3d

  Plot3D equivalent of face_center_from_volume.
  Plot3DZone stores x/y/z arrays as doubles with layout compatible with
  plot3d_index():
    idx = (k-1)*ni*nj + (j-1)*ni + (i-1)

  Parameters:
    - zone: Plot3DZone (idx/name, dims, x/y/z vectors)
    - face,u,v: same conventions as CGNS path
=====================================================================*/
Point face_center_from_plot3d(const Plot3DZone& zone, int face, int u, int v)
{
    bool is2D = (zone.nk() == 1);
    long long u0 = u + 1;
    long long v0 = v + 1;

    auto node = [&](long long i, long long j, long long k) {
        size_t idx = plot3d_index(zone, i, j, k);
        return Point{zone.x[idx], zone.y[idx], zone.z[idx]};
    };

    if (is2D) {
        if (face == 0 || face == 1) {
            long long j = u0;
            long long i = (face == 0) ? 1 : zone.ni();
            Point p1 = node(i, j, 1);
            Point p2 = node(i, j + 1, 1);
            return Point{(p1.x + p2.x) * 0.5, (p1.y + p2.y) * 0.5, (p1.z + p2.z) * 0.5};
        }
        if (face == 2 || face == 3) {
            long long i = u0;
            long long j = (face == 2) ? 1 : zone.nj();
            Point p1 = node(i, j, 1);
            Point p2 = node(i + 1, j, 1);
            return Point{(p1.x + p2.x) * 0.5, (p1.y + p2.y) * 0.5, (p1.z + p2.z) * 0.5};
        }
        return Point{0, 0, 0};
    }

    if (face == 0 || face == 1) {
        long long j = u0;
        long long k = v0;
        long long i = (face == 0) ? 1 : zone.ni();
        Point p1 = node(i, j, k);
        Point p2 = node(i, j + 1, k);
        Point p3 = node(i, j, k + 1);
        Point p4 = node(i, j + 1, k + 1);
        return Point{(p1.x + p2.x + p3.x + p4.x) * 0.25,
                     (p1.y + p2.y + p3.y + p4.y) * 0.25,
                     (p1.z + p2.z + p3.z + p4.z) * 0.25};
    }

    if (face == 2 || face == 3) {
        long long i = u0;
        long long k = v0;
        long long j = (face == 2) ? 1 : zone.nj();
        Point p1 = node(i, j, k);
        Point p2 = node(i + 1, j, k);
        Point p3 = node(i, j, k + 1);
        Point p4 = node(i + 1, j, k + 1);
        return Point{(p1.x + p2.x + p3.x + p4.x) * 0.25,
                     (p1.y + p2.y + p3.y + p4.y) * 0.25,
                     (p1.z + p2.z + p3.z + p4.z) * 0.25};
    }

    long long i = u0;
    long long j = v0;
    long long k = (face == 4) ? 1 : zone.nk();
    Point p1 = node(i, j, k);
    Point p2 = node(i + 1, j, k);
    Point p3 = node(i, j + 1, k);
    Point p4 = node(i + 1, j + 1, k);
    return Point{(p1.x + p2.x + p3.x + p4.x) * 0.25,
                 (p1.y + p2.y + p3.y + p4.y) * 0.25,
                 (p1.z + p2.z + p3.z + p4.z) * 0.25};
}

/*=====================================================================
  build_faces_cgns

  Main CGNS face extraction routine.

  For each Zone in mesh:
    - Decide whether to read full volume or per-face planes:
        use_volume = (vol_elems <= k * face_elems) && (vol_bytes <= MEM_CAP_BYTES)
      where k is K2D or K3D from k_values.hpp.

    - For each of the 6 face IDs:
        * compute face-local dimensions (uLen, vLen)
        * create a FaceGrid:
            - index:    int array of size uLen*vLen mapping (u,v)->FaceData index
            - consumed: uint8 array used later by resolve_matches to mark patches
        * compute center point for each face cell
        * append a FaceData entry via build.faces.push(...)

  FaceBuildResult layout notes:
    - build.faces is a Struct-of-Arrays:
        vectors of doubles x/y/z and ints block/face/u/v/i/j/k.
      push() returns the index into these arrays (int).
    - build.face_cells stores FaceData indices grouped by face_id.
      start/count per face_id are stored in build.cell_start/cell_count.

  Progress counting:
    - build.faces_loaded increments by 6 for 3D zones, 4 for 2D zones.
=====================================================================*/
FaceBuildResult build_faces_cgns(Mesh& mesh, bool show_progress)
{
    FaceBuildResult build;

    const int nBlocks = static_cast<int>(mesh.zones().size());

    // Storage is allocated for 6 faces per zone even if nk==1.
    // Unused faces in 2D are represented with uLen/vLen==0 and empty grids.
    size_t face_count = static_cast<size_t>(nBlocks) * 6;
    build.grids.resize(face_count);
    build.bounds.resize(face_count);
    build.cell_start.assign(face_count, 0);
    build.cell_count.assign(face_count, 0);
    build.face_sum.assign(face_count, Point{0, 0, 0});  // centroid accumulator
    build.face_counts.assign(face_count, 0);            // number of cells per face

    /*---------------------------------------------------------
      Estimate total number of face cells across all zones to
      reserve contiguous storage for build.faces and face_cells.
    ---------------------------------------------------------*/
    size_t total_cells_est = 0;
    for (const auto& z : mesh.zones()) {
        for (int face = 0; face < 6; ++face) {
            int uLen = 0, vLen = 0;
            face_dims(z, face, uLen, vLen);
            if (uLen > 0 && vLen > 0)
                total_cells_est += static_cast<size_t>(uLen) * vLen;
        }
    }
    build.face_cells.reserve(total_cells_est);
    build.faces.reserve(total_cells_est);

    const int fn = mesh.file_id();
    const int B  = mesh.base_id();

    for (int b = 0; b < nBlocks; ++b)
    {
        const Zone& z = mesh.zones()[b];

        /*-----------------------------------------------------
          Decide read strategy (plane-by-plane vs full volume).

          Computation uses doubles to avoid intermediate overflow
          when multiplying large ni/nj/nk values.
        -----------------------------------------------------*/
        const double k2d = K2D;
        const double k3d = K3D;

        const std::uint64_t mem_cap_bytes = MEM_CAP_BYTES;
        const double ni = static_cast<double>(z.ni());
        const double nj = static_cast<double>(z.nj());
        const double nk = static_cast<double>(z.nk());
        const double vol_elems = ni * nj * nk;

        const double face_elems = (z.nk() == 1)
            ? 2.0 * (ni + nj)
            : 2.0 * (ni * nj + ni * nk + nj * nk);

        const double k = (z.nk() == 1) ? k2d : k3d;

        const std::uint64_t vol_bytes = static_cast<std::uint64_t>(
            vol_elems * 3.0 * static_cast<double>(sizeof(double)));

        const bool use_volume = (vol_elems <= k * face_elems) &&
                                (vol_bytes <= mem_cap_bytes);

        VolumeCoords volume;
        if (use_volume) {
            volume = read_volume_coords(fn, B, z);
        }

        for (int face = 0; face < 6; ++face) {
            int uLen = 0, vLen = 0;
            face_dims(z, face, uLen, vLen);

            int face_id = b * 6 + face;
            FaceGrid& grid = build.grids[static_cast<size_t>(face_id)];

            // cell_start marks where this face's indices begin in build.face_cells.
            size_t start = build.face_cells.size();
            build.cell_start[static_cast<size_t>(face_id)] = start;

            /*-------------------------------------------------
              Faces with uLen/vLen==0 are skipped and left with
              empty index/consumed arrays. This happens for:
                - kMin/kMax in 2D zones
                - any degenerate zone dimension that yields <=1
                  vertex along a face direction.
            -------------------------------------------------*/
            if (uLen <= 0 || vLen <= 0) {
                grid.uLen = 0;
                grid.vLen = 0;
                grid.index.clear();
                grid.consumed.clear();
                build.cell_count[static_cast<size_t>(face_id)] = 0;
                continue;
            }

            // Plane-based read allocates face-local coordinate arrays per face.
            FacePlane plane;
            if (!use_volume) {
                plane = read_face_plane(fn, B, z, face);
            }

            grid.uLen = uLen;
            grid.vLen = vLen;

            // grid.index: int mapping (u,v) -> FaceData index (or -1).
            // Size = uLen*vLen.
            grid.index.assign(static_cast<size_t>(uLen) * vLen, -1);

            // grid.consumed: uint8 (0/1) for patch consumption during resolve_matches.
            grid.consumed.assign(static_cast<size_t>(uLen) * vLen, 0);

            /*-------------------------------------------------
              Populate one FaceData entry per face cell.
              FaceData indices are stable integers returned by
              build.faces.push(...).
            -------------------------------------------------*/
            for (int v = 0; v < vLen; ++v)
                for (int u = 0; u < uLen; ++u)
                {
                    // i/j/k are stored as 1-based CGNS indices.
                    int i = 0, j = 0, k = 0;
                    face_cell_indices(z, face, u, v, i, j, k);

                    Point center = use_volume
                        ? face_center_from_volume(volume, z, face, u, v)
                        : face_center_from_plane(plane, z, face, u, v);

                    expand_bounds(build.bounds[static_cast<size_t>(face_id)], center);
                    build.face_sum[static_cast<size_t>(face_id)].x += center.x;
                    build.face_sum[static_cast<size_t>(face_id)].y += center.y;
                    build.face_sum[static_cast<size_t>(face_id)].z += center.z;
                    build.face_counts[static_cast<size_t>(face_id)] += 1;

                    int idx = build.faces.push(b, face, u, v, i, j, k, center);

                    // Store the FaceData index in the (u,v) grid.
                    grid.index[face_grid_index(grid, u, v)] = idx;

                    // Store FaceData indices grouped by face.
                    build.face_cells.push_back(idx);
                }

            // cell_count is the number of FaceData entries appended for this face.
            build.cell_count[static_cast<size_t>(face_id)] =
                static_cast<int>(build.face_cells.size() - start);
        }

        if (show_progress) build.faces_loaded += (z.nk() > 1 ? 6 : 4);
    }

    return build;
}

/*=====================================================================
  build_faces_plot3d

  Plot3D equivalent of build_faces_cgns.
  No CGNS I/O is performed; face centers are derived directly from the
  Plot3DZone x/y/z arrays.

  Storage behavior matches the CGNS path:
    - 6 faces allocated per block
    - 2D blocks leave k faces empty (uLen/vLen==0)
=====================================================================*/
FaceBuildResult build_faces_plot3d(const std::vector<Plot3DZone>& zones)
{
    FaceBuildResult build;

    const int nBlocks = static_cast<int>(zones.size());
    size_t face_count = static_cast<size_t>(nBlocks) * 6;
    build.grids.resize(face_count);
    build.bounds.resize(face_count);
    build.cell_start.assign(face_count, 0);
    build.cell_count.assign(face_count, 0);
    build.face_sum.assign(face_count, Point{0, 0, 0});
    build.face_counts.assign(face_count, 0);

    size_t total_cells_est = 0;
    for (const auto& z : zones) {
        for (int face = 0; face < 6; ++face) {
            int uLen = 0, vLen = 0;
            face_dims(z, face, uLen, vLen);
            if (uLen > 0 && vLen > 0)
                total_cells_est += static_cast<size_t>(uLen) * vLen;
        }
    }
    build.face_cells.reserve(total_cells_est);
    build.faces.reserve(total_cells_est);

    for (int b = 0; b < nBlocks; ++b)
    {
        const Plot3DZone& z = zones[b];

        for (int face = 0; face < 6; ++face) {
            int uLen = 0, vLen = 0;
            face_dims(z, face, uLen, vLen);
            int face_id = b * 6 + face;
            FaceGrid& grid = build.grids[static_cast<size_t>(face_id)];
            size_t start = build.face_cells.size();
            build.cell_start[static_cast<size_t>(face_id)] = start;
            if (uLen <= 0 || vLen <= 0) {
                grid.uLen = 0;
                grid.vLen = 0;
                grid.index.clear();
                grid.consumed.clear();
                build.cell_count[static_cast<size_t>(face_id)] = 0;
                continue;
            }

            grid.uLen = uLen;
            grid.vLen = vLen;
            grid.index.assign(static_cast<size_t>(uLen) * vLen, -1);
            grid.consumed.assign(static_cast<size_t>(uLen) * vLen, 0);

            for (int v = 0; v < vLen; ++v)
                for (int u = 0; u < uLen; ++u)
                {
                    int i = 0, j = 0, k = 0;
                    face_cell_indices(z, face, u, v, i, j, k);
                    Point center = face_center_from_plot3d(z, face, u, v);
                    expand_bounds(build.bounds[static_cast<size_t>(face_id)], center);
                    build.face_sum[static_cast<size_t>(face_id)].x += center.x;
                    build.face_sum[static_cast<size_t>(face_id)].y += center.y;
                    build.face_sum[static_cast<size_t>(face_id)].z += center.z;
                    build.face_counts[static_cast<size_t>(face_id)] += 1;

                    int idx = build.faces.push(b, face, u, v, i, j, k, center);
                    grid.index[face_grid_index(grid, u, v)] = idx;
                    build.face_cells.push_back(idx);
                }

            build.cell_count[static_cast<size_t>(face_id)] =
                static_cast<int>(build.face_cells.size() - start);
        }

        build.faces_loaded += (z.nk() > 1 ? 6 : 4);
    }

    return build;
}

/*=====================================================================
  build_face_summaries

  Produces one FaceSummary per face_id (block*6 + face).

  FaceSummary contains:
    - identifiers (face_id, block, face)
    - uLen/vLen (face-local grid dimensions)
    - bounds (axis-aligned bounding box over all face-cell centers)
    - centroid (mean of face-cell centers)
    - normal (unit normal estimate) when uLen>1 and vLen>1

  Normal computation:
    - Uses the first 2x2 corner in face-local grid:
        p00 at (0,0)
        p10 at (1,0)
        p01 at (0,1)
    - v1 = p10 - p00, v2 = p01 - p00
    - n = cross(v1, v2), normalized if non-zero
  This is a local orientation estimate used in later matching logic.
=====================================================================*/
std::vector<FaceSummary> build_face_summaries(const FaceBuildResult& build)
{
    size_t face_count = build.grids.size();
    std::vector<FaceSummary> summaries(face_count);

    // Pass 1: populate basic metadata + centroid from accumulated sums.
    for (size_t face_id = 0; face_id < face_count; ++face_id) {
        FaceSummary summary;
        summary.face_id = static_cast<int>(face_id);
        summary.block = static_cast<int>(face_id / 6);
        summary.face = static_cast<int>(face_id % 6);
        summary.uLen = build.grids[face_id].uLen;
        summary.vLen = build.grids[face_id].vLen;
        summary.bounds = build.bounds[face_id];
        if (build.face_counts[face_id] > 0) {
            summary.centroid = Point{
                build.face_sum[face_id].x / build.face_counts[face_id],
                build.face_sum[face_id].y / build.face_counts[face_id],
                build.face_sum[face_id].z / build.face_counts[face_id]
            };
        }
        summaries[face_id] = summary;
    }

    // Pass 2: compute an estimated face normal when the face has at least a 2x2 grid.
    for (size_t face_id = 0; face_id < face_count; ++face_id) {
        FaceSummary& summary = summaries[face_id];
        if (summary.uLen > 1 && summary.vLen > 1) {
            const FaceGrid& grid = build.grids[face_id];

            // FaceGrid index holds FaceData indices, which index into build.faces SoA arrays.
            int idx00 = grid.index[face_grid_index(grid, 0, 0)];
            int idx10 = grid.index[face_grid_index(grid, 1, 0)];
            int idx01 = grid.index[face_grid_index(grid, 0, 1)];

            Point p00{build.faces.x[idx00], build.faces.y[idx00], build.faces.z[idx00]};
            Point p10{build.faces.x[idx10], build.faces.y[idx10], build.faces.z[idx10]};
            Point p01{build.faces.x[idx01], build.faces.y[idx01], build.faces.z[idx01]};

            Point v1{p10.x - p00.x, p10.y - p00.y, p10.z - p00.z};
            Point v2{p01.x - p00.x, p01.y - p00.y, p01.z - p00.z};
            
            Point n = cross(v1, v2);
            double nlen = norm(n);
            if (nlen > 0.0) {
                summary.normal = normalize(n);
                summary.normal_valid = true;
            }
        }
    }

    return summaries;
}

} // namespace fs::conn
