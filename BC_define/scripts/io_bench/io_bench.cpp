// Minimal CGNS IO benchmark tool for k calibration.
#include <cgnslib.h>

#include <chrono>
#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>

#ifndef CGNS_MAX_NAME_LENGTH
#define CGNS_MAX_NAME_LENGTH 32
#endif

static void usage()
{
    std::cerr << "Usage: io_bench mesh.cgns [--bench-iter N]\n";
}

int main(int argc, char** argv)
{
    if (argc < 2) {
        usage();
        return 1;
    }

    std::string mesh = argv[1];
    int iters = 3;
    for (int i = 2; i < argc; ++i) {
        std::string flag = argv[i];
        if (flag == "--bench-iter") {
            if (i + 1 >= argc) {
                usage();
                return 1;
            }
            iters = std::atoi(argv[++i]);
        }
    }

    int fn = 0;
    if (cg_open(mesh.c_str(), CG_MODE_READ, &fn)) {
        cg_error_print();
        return 1;
    }

    int B = 1;
    int nZones = 0;
    if (cg_nzones(fn, B, &nZones)) {
        cg_error_print();
        cg_close(fn);
        return 1;
    }

    static constexpr const char* cname[3] =
        {"CoordinateX", "CoordinateY", "CoordinateZ"};

    std::cout << "IO benchmark (iters=" << iters << ")\n";

    double total_vol_time = 0.0;
    double total_face_time = 0.0;
    long long total_vol_elems = 0;
    long long total_face_elems = 0;

    for (int Z = 1; Z <= nZones; ++Z) {
        char zname[CGNS_MAX_NAME_LENGTH + 1] = "";
        cgsize_t size[3][3];
        if (cg_zone_read(fn, B, Z, zname, &size[0][0])) {
            cg_error_print();
            cg_close(fn);
            return 1;
        }

        long long ni = size[0][0];
        long long nj = size[0][1];
        long long nk = size[0][2];
        if (nk < 1) nk = 1;

        long long vol_elems = ni * nj * nk;
        if (vol_elems <= 0)
            continue;

        struct FaceRead {
            cgsize_t s[3]{};
            cgsize_t e[3]{};
            long long count = 0;
        };

        FaceRead faces[6];
        long long face_elems = 0;
        long long max_face = 0;

        for (int face = 0; face < 6; ++face) {
            cgsize_t s[3] = {1, 1, 1};
            cgsize_t e[3] = {static_cast<cgsize_t>(ni),
                             static_cast<cgsize_t>(nj),
                             static_cast<cgsize_t>(nk)};
            long long dim1 = 0;
            long long dim2 = 0;

            if (face == 0) {
                s[0] = 1; e[0] = 1; dim1 = nj; dim2 = nk;
            } else if (face == 1) {
                s[0] = static_cast<cgsize_t>(ni);
                e[0] = static_cast<cgsize_t>(ni);
                dim1 = nj; dim2 = nk;
            } else if (face == 2) {
                s[1] = 1; e[1] = 1; dim1 = ni; dim2 = nk;
            } else if (face == 3) {
                s[1] = static_cast<cgsize_t>(nj);
                e[1] = static_cast<cgsize_t>(nj);
                dim1 = ni; dim2 = nk;
            } else if (face == 4) {
                s[2] = 1; e[2] = 1; dim1 = ni; dim2 = nj;
            } else {
                s[2] = static_cast<cgsize_t>(nk);
                e[2] = static_cast<cgsize_t>(nk);
                dim1 = ni; dim2 = nj;
            }

            long long count = dim1 * dim2;
            if (count <= 0) {
                faces[face].count = 0;
                continue;
            }
            faces[face] = FaceRead{{s[0], s[1], s[2]}, {e[0], e[1], e[2]}, count};
            face_elems += count;
            if (count > max_face)
                max_face = count;
        }

        if (face_elems <= 0)
            continue;

        std::vector<double> vol_buf(static_cast<size_t>(vol_elems));
        std::vector<double> face_buf(static_cast<size_t>(max_face));
        cgsize_t sv[3] = {1, 1, 1};
        cgsize_t ev[3] = {static_cast<cgsize_t>(ni),
                          static_cast<cgsize_t>(nj),
                          static_cast<cgsize_t>(nk)};

        auto read_volume = [&]() {
            for (int c = 0; c < 3; ++c) {
                if (cg_coord_read(fn, B, Z, cname[c],
                                  CGNS_ENUMV(RealDouble), sv, ev,
                                  vol_buf.data())) {
                    cg_error_print();
                    return false;
                }
            }
            return true;
        };

        auto read_faces = [&]() {
            for (int face = 0; face < 6; ++face) {
                if (faces[face].count <= 0)
                    continue;
                for (int c = 0; c < 3; ++c) {
                    if (cg_coord_read(fn, B, Z, cname[c],
                                      CGNS_ENUMV(RealDouble),
                                      faces[face].s, faces[face].e,
                                      face_buf.data())) {
                        cg_error_print();
                        return false;
                    }
                }
            }
            return true;
        };

        if (!read_volume() || !read_faces()) {
            cg_close(fn);
            return 1;
        }

        auto t0 = std::chrono::steady_clock::now();
        for (int i = 0; i < iters; ++i)
            if (!read_volume()) return 1;
        auto t1 = std::chrono::steady_clock::now();
        for (int i = 0; i < iters; ++i)
            if (!read_faces()) return 1;
        auto t2 = std::chrono::steady_clock::now();

        std::chrono::duration<double> vol_dur = t1 - t0;
        std::chrono::duration<double> face_dur = t2 - t1;

        double vol_time = vol_dur.count() / iters;
        double face_time = face_dur.count() / iters;
        double t_v = vol_time / static_cast<double>(vol_elems);
        double t_f = face_time / static_cast<double>(face_elems);
        double k_face = (t_v > 0.0) ? (t_f / t_v) : 0.0;

        std::cout << "Zone " << Z << " (" << zname << ") "
                  << ni << "x" << nj << "x" << nk << "\n"
                  << "  vol elems=" << vol_elems
                  << " face elems=" << face_elems
                  << " Tvol=" << (vol_time * 1e3) << " ms"
                  << " Tface=" << (face_time * 1e3) << " ms\n"
                  << "  t_v=" << (t_v * 1e9) << " ns/elem"
                  << " t_f=" << (t_f * 1e9) << " ns/elem"
                  << " k_face=" << k_face << "\n";

        total_vol_time += vol_time;
        total_face_time += face_time;
        total_vol_elems += vol_elems;
        total_face_elems += face_elems;
    }

    if (total_vol_elems > 0 && total_face_elems > 0) {
        double t_v = total_vol_time / static_cast<double>(total_vol_elems);
        double t_f = total_face_time / static_cast<double>(total_face_elems);
        double k_face = (t_v > 0.0) ? (t_f / t_v) : 0.0;
        std::cout << "Overall: t_v=" << (t_v * 1e9) << " ns/elem"
                  << " t_f=" << (t_f * 1e9) << " ns/elem"
                  << " k_face=" << k_face << "\n"
                  << "Rule: use volume if V <= k_face * (face_elems)\n";
    }

    cg_close(fn);
    return 0;
}
