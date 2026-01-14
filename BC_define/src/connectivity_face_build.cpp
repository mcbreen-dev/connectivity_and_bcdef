//─────────────────────────────────────────────────────────────
// File: src/connectivity_face_build.cpp
// Face extraction and metadata construction
//─────────────────────────────────────────────────────────────
#include "connectivity_internal.hpp"
#include "k_values.hpp"

#include <algorithm>

extern "C" {
    #include <cgnslib.h>
}

namespace fs::conn {

Point face_center_from_plane(const FacePlane& plane, const Zone& z, int face, int u, int v)
{
    bool is2D = (z.nk() == 1);
    long long u0 = u + 1;
    long long v0 = v + 1;

    if (is2D) {
        if (face == 0 || face == 1) {
            long long j = u0;
            size_t p1 = plane_index(plane, j, 1);
            size_t p2 = plane_index(plane, j + 1, 1);
            return Point{(plane.x[p1] + plane.x[p2]) * 0.5,
                         (plane.y[p1] + plane.y[p2]) * 0.5,
                         (plane.z[p1] + plane.z[p2]) * 0.5};
        }
        if (face == 2 || face == 3) {
            long long i = u0;
            size_t p1 = plane_index(plane, i, 1);
            size_t p2 = plane_index(plane, i + 1, 1);
            return Point{(plane.x[p1] + plane.x[p2]) * 0.5,
                         (plane.y[p1] + plane.y[p2]) * 0.5,
                         (plane.z[p1] + plane.z[p2]) * 0.5};
        }
        return Point{0, 0, 0};
    }

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

FacePlane read_face_plane(int fn, int B, const Zone& z, int face)
{
    FacePlane plane;
    plane.face = face;

    cgsize_t s[3] = {1, 1, 1};
    cgsize_t e[3] = {static_cast<cgsize_t>(z.ni()),
                     static_cast<cgsize_t>(z.nj()),
                     static_cast<cgsize_t>(z.nk())};

    if (face == 0) { // iMin
        s[0] = 1;
        e[0] = 1;
        plane.dim1 = z.nj();
        plane.dim2 = z.nk();
    } else if (face == 1) { // iMax
        s[0] = static_cast<cgsize_t>(z.ni());
        e[0] = static_cast<cgsize_t>(z.ni());
        plane.dim1 = z.nj();
        plane.dim2 = z.nk();
    } else if (face == 2) { // jMin
        s[1] = 1;
        e[1] = 1;
        plane.dim1 = z.ni();
        plane.dim2 = z.nk();
    } else if (face == 3) { // jMax
        s[1] = static_cast<cgsize_t>(z.nj());
        e[1] = static_cast<cgsize_t>(z.nj());
        plane.dim1 = z.ni();
        plane.dim2 = z.nk();
    } else if (face == 4) { // kMin
        s[2] = 1;
        e[2] = 1;
        plane.dim1 = z.ni();
        plane.dim2 = z.nj();
    } else { // kMax
        s[2] = static_cast<cgsize_t>(z.nk());
        e[2] = static_cast<cgsize_t>(z.nk());
        plane.dim1 = z.ni();
        plane.dim2 = z.nj();
    }

    if (plane.dim1 <= 0 || plane.dim2 <= 0)
        return plane;

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

FaceBuildResult build_faces_cgns(Mesh& mesh, bool show_progress)
{
    FaceBuildResult build;

    const int nBlocks = static_cast<int>(mesh.zones().size());
    size_t face_count = static_cast<size_t>(nBlocks) * 6;
    build.grids.resize(face_count);
    build.bounds.resize(face_count);
    build.cell_start.assign(face_count, 0);
    build.cell_count.assign(face_count, 0);
    build.face_sum.assign(face_count, Point{0, 0, 0});
    build.face_counts.assign(face_count, 0);

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

            FacePlane plane;
            if (!use_volume) {
                plane = read_face_plane(fn, B, z, face);
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
                    Point center = use_volume
                        ? face_center_from_volume(volume, z, face, u, v)
                        : face_center_from_plane(plane, z, face, u, v);
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

        if (show_progress) build.faces_loaded += (z.nk() > 1 ? 6 : 4);
    }

    return build;
}

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

std::vector<FaceSummary> build_face_summaries(const FaceBuildResult& build)
{
    size_t face_count = build.grids.size();
    std::vector<FaceSummary> summaries(face_count);

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

    for (size_t face_id = 0; face_id < face_count; ++face_id) {
        FaceSummary& summary = summaries[face_id];
        if (summary.uLen > 1 && summary.vLen > 1) {
            const FaceGrid& grid = build.grids[face_id];
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
