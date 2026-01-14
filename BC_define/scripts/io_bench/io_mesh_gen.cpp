// Simple structured CGNS mesh generator for IO benchmarking.
#include <cgnslib.h>

#include <array>
#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>

struct AxisSplit {
    std::vector<int> start;
    std::vector<int> count;
};

static AxisSplit split_axis(int n, int parts)
{
    AxisSplit out;
    if (parts < 1 || n < 1) return out;

    int total_cells = n - 1;
    int base = total_cells / parts;
    int rem = total_cells % parts;
    int start_cell = 0;

    out.start.reserve(parts);
    out.count.reserve(parts);

    for (int p = 0; p < parts; ++p) {
        int cells = base + (p < rem ? 1 : 0);
        int verts = cells + 1;
        int start_vertex = start_cell + 1;
        out.start.push_back(start_vertex);
        out.count.push_back(verts);
        start_cell += cells;
    }

    return out;
}

static void usage()
{
    std::cerr << "Usage: io_mesh_gen out.cgns ni nj nk [sx sy sz]\n";
}

int main(int argc, char** argv)
{
    if (argc < 5) {
        usage();
        return 1;
    }

    const std::string out_path = argv[1];
    const int ni = std::atoi(argv[2]);
    const int nj = std::atoi(argv[3]);
    const int nk = std::atoi(argv[4]);
    int sx = 1, sy = 1, sz = 1;
    if (argc >= 8) {
        sx = std::atoi(argv[5]);
        sy = std::atoi(argv[6]);
        sz = std::atoi(argv[7]);
    }

    if (ni < 2 || nj < 2 || nk < 1 || sx < 1 || sy < 1 || sz < 1) {
        std::cerr << "Invalid dimensions or splits\n";
        return 1;
    }

    AxisSplit xs = split_axis(ni, sx);
    AxisSplit ys = split_axis(nj, sy);
    AxisSplit zs = split_axis(nk, sz);

    int fn = 0;
    if (cg_open(out_path.c_str(), CG_MODE_WRITE, &fn)) {
        cg_error_print();
        return 1;
    }

    int B = 0;
    if (cg_base_write(fn, "Base", 3, 3, &B)) {
        cg_error_print();
        cg_close(fn);
        return 1;
    }

    int zone_id = 0;
    for (int kz = 0; kz < sz; ++kz) {
        for (int jy = 0; jy < sy; ++jy) {
            for (int ix = 0; ix < sx; ++ix) {
                ++zone_id;
                int ni_b = xs.count[ix];
                int nj_b = ys.count[jy];
                int nk_b = zs.count[kz];

                cgsize_t size[9];
                size[0] = ni_b;
                size[1] = nj_b;
                size[2] = nk_b;
                size[3] = ni_b - 1;
                size[4] = nj_b - 1;
                size[5] = nk_b - 1;
                size[6] = ni_b;
                size[7] = nj_b;
                size[8] = nk_b;

                std::string zname = "block-" + std::to_string(zone_id);
                int Z = 0;
                if (cg_zone_write(fn, B, zname.c_str(), size,
                                  CGNS_ENUMV(Structured), &Z)) {
                    cg_error_print();
                    cg_close(fn);
                    return 1;
                }

                int i0 = xs.start[ix] - 1;
                int j0 = ys.start[jy] - 1;
                int k0 = zs.start[kz] - 1;

                std::vector<double> xplane(static_cast<size_t>(ni_b) * nj_b);
                std::vector<double> yplane(static_cast<size_t>(ni_b) * nj_b);
                std::vector<double> zplane(static_cast<size_t>(ni_b) * nj_b);

                for (int k = 0; k < nk_b; ++k) {
                    size_t idx = 0;
                    for (int j = 0; j < nj_b; ++j) {
                        double y = static_cast<double>(j0 + j);
                        for (int i = 0; i < ni_b; ++i) {
                            xplane[idx] = static_cast<double>(i0 + i);
                            yplane[idx] = y;
                            zplane[idx] = static_cast<double>(k0 + k);
                            ++idx;
                        }
                    }

                    cgsize_t s[3] = {1, 1, static_cast<cgsize_t>(k + 1)};
                    cgsize_t e[3] = {static_cast<cgsize_t>(ni_b),
                                     static_cast<cgsize_t>(nj_b),
                                     static_cast<cgsize_t>(k + 1)};
                    int C = 0;
                    if (cg_coord_partial_write(fn, B, Z, CGNS_ENUMV(RealDouble),
                                               "CoordinateX",
                                               s, e, xplane.data(), &C) ||
                        cg_coord_partial_write(fn, B, Z, CGNS_ENUMV(RealDouble),
                                               "CoordinateY",
                                               s, e, yplane.data(), &C) ||
                        cg_coord_partial_write(fn, B, Z, CGNS_ENUMV(RealDouble),
                                               "CoordinateZ",
                                               s, e, zplane.data(), &C)) {
                        cg_error_print();
                        cg_close(fn);
                        return 1;
                    }
                }
            }
        }
    }

    cg_close(fn);
    return 0;
}
