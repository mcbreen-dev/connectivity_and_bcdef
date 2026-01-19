//─────────────────────────────────────────────────────────────
// File: src/bc_define_cgns.cpp
// CGNS mesh I/O and boundary-condition writing
//─────────────────────────────────────────────────────────────
//
// This file contains the CGNS-facing portion of the BC_define pipeline.
// It is responsible for translating resolved boundary patches and
// user-defined BC specifications into CGNS BC_t nodes written
// into the mesh file.
//
// Major responsibilities:
//   - Benchmark CGNS coordinate I/O performance (optional diagnostic).
//   - Map face directions to CGNS GridLocation_t values.
//   - Detect which structured face a CGNS PointRange lies on.
//   - Ensure CGNS Family_t nodes exist and carry the correct BC type.
//   - Safely write BC_t nodes while avoiding duplicates and overlaps.
//   - Enforce strict validation for user-specified BCs while allowing
//     conservative behavior for automatically generated BCs.
//
// Design philosophy:
//   - Be conservative with existing CGNS data: never overwrite or
//     partially modify BCs in place.
//   - Detect and error out on conflicting user specifications.
//   - Allow auto-generated BCs (autowall / autofarfield) to yield to
//     any existing BC definitions.
//   - Cache CGNS lookups aggressively to minimize repeated I/O.
//
// This file does *not* perform connectivity detection or BC parsing;
// those tasks are handled by connectivity_* and bc_define_bcspec.cpp
// respectively.
//─────────────────────────────────────────────────────────────
#include "bc_define_internal.hpp"
#include "logger.hpp"
#include "cgns_cleanup.hpp"

#include <chrono>
#include <iostream>
#include <unordered_map>

namespace bcdef::boundary {

/*=====================================================================
  face_dims (local helper)

  Compute the face-local cell grid dimensions (uLen, vLen) for a given
  structured zone face.

  This is a simplified, BC-local equivalent of the face dimension logic
  used in the connectivity pipeline.

  2-D behavior:
    - For Nk == 1, only i- and j-faces are valid.
    - vLen is forced to 1 so the rest of the logic can treat faces as
      degenerate 2-D grids.
=====================================================================*/
static void face_dims(const bcdef::Zone& z, int face, int& uLen, int& vLen)
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
  run_io_benchmark

  Optional diagnostic routine that measures CGNS coordinate read
  performance for:
    - full-volume reads
    - per-face reads

  The benchmark reports:
    - time per element for volume vs face reads
    - relative cost factor (k_face)
    - a heuristic indicating whether volume or face reads are cheaper

  This is intended for developer diagnostics only and has no effect
  on BC writing logic.
=====================================================================*/
void run_io_benchmark(bcdef::Mesh& mesh, bcdef::Logger& log, int iters)
{
    const int fn = mesh.file_id();
    const int B  = mesh.base_id();
    static constexpr const char* cname[3] =
        {"CoordinateX", "CoordinateY", "CoordinateZ"};

    log.info("IO benchmark started");
    std::cout << "IO benchmark (iters=" << iters << ")\n";

    double total_vol_time = 0.0;
    double total_face_time = 0.0;
    long long total_vol_elems = 0;
    long long total_face_elems = 0;

    for (const auto& z : mesh.zones()) {
        long long ni = z.ni();
        long long nj = z.nj();
        long long nk = z.nk();
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

        /*---------------------------------------------------------
          Prepare per-face coordinate ranges
        ---------------------------------------------------------*/
        for (int face = 0; face < 6; ++face) {
            int uLen = 0, vLen = 0;
            face_dims(z, face, uLen, vLen);
            if (uLen <= 0 || vLen <= 0) {
                faces[face].count = 0;
                continue;
            }

            cgsize_t s[3] = {1, 1, 1};
            cgsize_t e[3] = {static_cast<cgsize_t>(ni),
                             static_cast<cgsize_t>(nj),
                             static_cast<cgsize_t>(nk)};
            long long dim1 = 0;
            long long dim2 = 0;

            if (face == 0) { // iMin
                s[0] = 1;
                e[0] = 1;
                dim1 = nj;
                dim2 = nk;
            } else if (face == 1) { // iMax
                s[0] = static_cast<cgsize_t>(ni);
                e[0] = static_cast<cgsize_t>(ni);
                dim1 = nj;
                dim2 = nk;
            } else if (face == 2) { // jMin
                s[1] = 1;
                e[1] = 1;
                dim1 = ni;
                dim2 = nk;
            } else if (face == 3) { // jMax
                s[1] = static_cast<cgsize_t>(nj);
                e[1] = static_cast<cgsize_t>(nj);
                dim1 = ni;
                dim2 = nk;
            } else if (face == 4) { // kMin
                s[2] = 1;
                e[2] = 1;
                dim1 = ni;
                dim2 = nj;
            } else { // kMax
                s[2] = static_cast<cgsize_t>(nk);
                e[2] = static_cast<cgsize_t>(nk);
                dim1 = ni;
                dim2 = nj;
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
                CG_CALL(cg_coord_read(fn, B, z.idx, cname[c],
                                      CGNS_ENUMV(RealDouble), sv, ev,
                                      vol_buf.data()),
                        "cg_coord_read volume failed");
            }
        };

        auto read_faces = [&]() {
            for (int face = 0; face < 6; ++face) {
                if (faces[face].count <= 0)
                    continue;
                for (int c = 0; c < 3; ++c) {
                    CG_CALL(cg_coord_read(fn, B, z.idx, cname[c],
                                          CGNS_ENUMV(RealDouble),
                                          faces[face].s, faces[face].e,
                                          face_buf.data()),
                            "cg_coord_read face failed");
                }
            }
        };

        read_volume();
        read_faces();

        auto t0 = std::chrono::steady_clock::now();
        for (int i = 0; i < iters; ++i)
            read_volume();
        auto t1 = std::chrono::steady_clock::now();
        for (int i = 0; i < iters; ++i)
            read_faces();
        auto t2 = std::chrono::steady_clock::now();

        std::chrono::duration<double> vol_dur = t1 - t0;
        std::chrono::duration<double> face_dur = t2 - t1;

        double vol_time = vol_dur.count() / iters;
        double face_time = face_dur.count() / iters;
        double t_v = vol_time / static_cast<double>(vol_elems);
        double t_f = face_time / static_cast<double>(face_elems);
        double k_face = (t_v > 0.0) ? (t_f / t_v) : 0.0;

        long long surf_base = (nk == 1) ? (ni + nj)
                                        : (ni * nj + ni * nk + nj * nk);
        double k_surface = (surf_base > 0)
            ? (k_face * static_cast<double>(face_elems) / static_cast<double>(surf_base))
            : 0.0;

        std::cout << "Zone " << z.idx << " (" << z.name << ") "
                  << ni << "x" << nj << "x" << nk << "\n"
                  << "  vol elems=" << vol_elems
                  << " face elems=" << face_elems
                  << " Tvol=" << (vol_time * 1e3) << " ms"
                  << " Tface=" << (face_time * 1e3) << " ms\n"
                  << "  t_v=" << (t_v * 1e9) << " ns/elem"
                  << " t_f=" << (t_f * 1e9) << " ns/elem"
                  << " k_face=" << k_face
                  << " k_surface=" << k_surface << "\n";

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

    log.info("IO benchmark finished");
}

/*=====================================================================
  face_grid_location

  Map a structured face direction to the appropriate CGNS GridLocation_t.

  Behavior:
    - 3-D meshes: returns IFaceCenter / JFaceCenter / KFaceCenter
    - 2-D meshes: returns EdgeCenter
=====================================================================*/
CGNS_ENUMT(GridLocation_t) face_grid_location(bcdef::FaceDir face, int cell_dim)
{
    if (cell_dim >= 3) {
        switch (face) {
            case bcdef::FaceDir::IMIN:
            case bcdef::FaceDir::IMAX:
                return CGNS_ENUMV(IFaceCenter);
            case bcdef::FaceDir::JMIN:
            case bcdef::FaceDir::JMAX:
                return CGNS_ENUMV(JFaceCenter);
            case bcdef::FaceDir::KMIN:
            case bcdef::FaceDir::KMAX:
                return CGNS_ENUMV(KFaceCenter);
        }
        return CGNS_ENUMV(FaceCenter);
    }
    if (cell_dim == 2)
        return CGNS_ENUMV(EdgeCenter);
    throw std::runtime_error("Unsupported CellDimension for BC grid location");
}

/*=====================================================================
  range_face

  Determine whether a CGNS PointRange corresponds to a full structured
  face of the given zone, and if so, which face.

  Returns:
    - true  and sets 'out' if the range lies on a face
    - false otherwise
=====================================================================*/
bool range_face(const bcdef::Zone& zone,
                const std::array<cgsize_t, 6>& r,
                bcdef::FaceDir& out)
{
    if (r[0] == r[3]) {
        if (r[0] == 1) {
            out = bcdef::FaceDir::IMIN;
            return true;
        }
        if (r[0] == static_cast<cgsize_t>(zone.ni() - 1)) {
            out = bcdef::FaceDir::IMAX;
            return true;
        }
        return false;
    }
    if (r[1] == r[4]) {
        if (r[1] == 1) {
            out = bcdef::FaceDir::JMIN;
            return true;
        }
        if (r[1] == static_cast<cgsize_t>(zone.nj() - 1)) {
            out = bcdef::FaceDir::JMAX;
            return true;
        }
        return false;
    }
    if (zone.nk() == 1)
        return false;
    if (r[2] == r[5]) {
        if (r[2] == 1) {
            out = bcdef::FaceDir::KMIN;
            return true;
        }
        if (r[2] == static_cast<cgsize_t>(zone.nk() - 1)) {
            out = bcdef::FaceDir::KMAX;
            return true;
        }
        return false;
    }
    return false;
}

/*=====================================================================
  ensure_family

  Ensure that a CGNS Family_t with the given name exists and has an
  associated FamilyBC of the requested BC type.

  Uses a cache to avoid repeated CGNS tree traversal.
=====================================================================*/
int ensure_family(int fn, int B, const std::string& family,
                  CGNS_ENUMT(BCType_t) bc_type,
                  std::unordered_map<std::string, int>& cache)
{
    auto it = cache.find(family);
    if (it != cache.end())
        return it->second;

    int nFamilies = 0;
    CG_CALL(cg_nfamilies(fn, B, &nFamilies), "cg_nfamilies failed");
    char name[CGNS_MAX_NAME_LENGTH + 1];
    for (int f = 1; f <= nFamilies; ++f) {
        int nboco = 0;
        int ngeos = 0;
        CG_CALL(cg_family_read(fn, B, f, name, &nboco, &ngeos),
                "cg_family_read failed");
        if (family == name) {
            if (nboco == 0) {
                int bc = 0;
                CG_CALL(cg_fambc_write(fn, B, f, "FamilyBC", bc_type, &bc),
                        "cg_fambc_write failed");
            }
            cache[family] = f;
            return f;
        }
    }

    int F = 0;
    CG_CALL(cg_family_write(fn, B, family.c_str(), &F),
            "cg_family_write failed");
    int bc = 0;
    CG_CALL(cg_fambc_write(fn, B, F, "FamilyBC", bc_type, &bc),
            "cg_fambc_write failed");
    cache[family] = F;
    return F;
}

/*=====================================================================
  write_boundary_conditions

  Main entry point for BC emission.

  Inputs:
    - mesh        : CGNS mesh (open for modification)
    - patches     : boundary patches detected by connectivity
    - specs       : user-defined BC specifications (per zone/face)
    - autowall    : assign wall BCs to unspecified faces
    - autofarfield: assign farfield BCs to unspecified faces

  Behavior:
    - Validates against existing BCs to prevent conflicts.
    - Writes new BC_t nodes only when safe.
    - Uses FamilySpecified BCs with FamilyName linkage.
=====================================================================*/
void write_boundary_conditions(
    bcdef::Mesh& mesh,
    const std::vector<bcdef::BoundaryPatch>& patches,
    const std::unordered_map<ZoneFaceKey, BCSpec, ZoneFaceKeyHash>& specs,
    bool autowall,
    bool autofarfield)
{
    int fn = mesh.file_id();
    int B = mesh.base_id();

    std::unordered_map<std::string, int> family_cache;
    std::unordered_map<ZoneFaceKey, int, ZoneFaceKeyHash> patch_counts;

    struct ExistingBC {
        std::string name;
        std::string family;
        CGNS_ENUMT(BCType_t) type = CGNS_ENUMV(BCTypeUserDefined);
        std::array<cgsize_t, 6> range{};
        bool has_range = false;
    };
    struct ZoneBCIndex {
        std::vector<ExistingBC> entries;
        std::unordered_map<std::string, size_t> by_name;
    };

    std::unordered_map<int, ZoneBCIndex> existing_by_zone;

    auto to_range = [&](const bcdef::PointRange& pr) {
        std::array<cgsize_t, 6> out{};
        if (mesh.cell_dim() == 2) {
            out[0] = static_cast<cgsize_t>(pr.begin[0]);
            out[1] = static_cast<cgsize_t>(pr.begin[1]);
            out[2] = 1;
            out[3] = static_cast<cgsize_t>(pr.end[0]);
            out[4] = static_cast<cgsize_t>(pr.end[1]);
            out[5] = 1;
        } else {
            out[0] = static_cast<cgsize_t>(pr.begin[0]);
            out[1] = static_cast<cgsize_t>(pr.begin[1]);
            out[2] = static_cast<cgsize_t>(pr.begin[2]);
            out[3] = static_cast<cgsize_t>(pr.end[0]);
            out[4] = static_cast<cgsize_t>(pr.end[1]);
            out[5] = static_cast<cgsize_t>(pr.end[2]);
        }
        return out;
    };

    for (const auto& zone : mesh.zones()) {
        int Z = zone.idx;
        int nBC = 0;
        if (cg_nbocos(fn, B, Z, &nBC) != CG_OK)
            continue;

        ZoneBCIndex index;
        for (int ibc = 1; ibc <= nBC; ++ibc) {
            char name[CGNS_MAX_NAME_LENGTH + 1] = {};
            CGNS_ENUMT(BCType_t) type;
            CGNS_ENUMT(PointSetType_t) pset;
            cgsize_t nPnts = 0;
            int NormalInd[3] = {0, 0, 0};
            cgsize_t NormalListSz = 0;
            CGNS_ENUMT(DataType_t) NormalDataType;
            int nDataset = 0;

            if (cg_boco_info(fn, B, Z, ibc, name, &type, &pset, &nPnts,
                             NormalInd, &NormalListSz, &NormalDataType, &nDataset) != CG_OK)
                continue;

            ExistingBC entry;
            entry.name = name;
            entry.type = type;
            if (type == CGNS_ENUMV(FamilySpecified)) {
                if (cg_goto(fn, B, "Zone_t", Z, "ZoneBC_t", 1,
                            "BC_t", ibc, "end") == CG_OK) {
                    char fam[CGNS_MAX_NAME_LENGTH + 1] = {};
                    if (cg_famname_read(fam) == CG_OK)
                        entry.family = fam;
                }
            }
            if (pset == CGNS_ENUMV(::PointRange) && nPnts == 2) {
                cgsize_t rng[6] = {0};
                if (cg_boco_read(fn, B, Z, ibc, rng, nullptr) == CG_OK) {
                    if (mesh.cell_dim() == 2) {
                        entry.range = {rng[0], rng[1], 1, rng[2], rng[3], 1};
                    } else {
                        entry.range = {rng[0], rng[1], rng[2], rng[3], rng[4], rng[5]};
                    }
                    entry.has_range = true;
                }
            }

            index.by_name[entry.name] = index.entries.size();
            index.entries.push_back(entry);
        }

        existing_by_zone[Z] = std::move(index);
    }

    for (const auto& patch : patches) {
        ZoneFaceKey key{patch.zone, patch.face};
        auto it = specs.find(key);
        const BCSpec* spec = nullptr;
        BCSpec auto_spec{};
        auto& zone_index = existing_by_zone[patch.zone];

        if (it == specs.end()) {
            if (!(autowall || autofarfield))
                continue;
            if (autowall) {
                auto_spec.family = "wall";
                auto_spec.type = CGNS_ENUMV(BCWall);
            } else {
                auto_spec.family = "farfield";
                auto_spec.type = CGNS_ENUMV(BCFarfield);
            }
            spec = &auto_spec;
        } else {
            spec = &it->second;
        }

        const bcdef::Zone& zone = mesh.zones()[patch.zone - 1];
        bcdef::PointRange face_range = face_center_range(zone, patch.face,
                                                      patch.vtxBegin, patch.vtxEnd);
        std::array<cgsize_t, 6> new_range = to_range(face_range);
        bool user_spec = (it != specs.end());

        bool overlap_skip = false;
        for (const auto& existing : zone_index.entries) {
            if (!existing.has_range)
                continue;
            bcdef::FaceDir existing_face;
            if (!range_face(zone, existing.range, existing_face))
                continue;
            if (existing_face != patch.face)
                continue;
            std::array<cgsize_t, 6> overlap{};
            if (!range_overlap(existing.range, new_range, overlap))
                continue;

            bool same_bc = (!existing.family.empty() &&
                            existing.family == spec->family &&
                            existing.type == spec->type);
            bool exact = ranges_equal(existing.range, new_range);

            if (user_spec) {
                if (same_bc && exact) {
                    overlap_skip = true;
                    break;
                }
                std::string existing_desc = existing.name;
                if (!existing.family.empty())
                    existing_desc = existing.family;
                throw std::runtime_error(
                    "Overlapping BCs in zone " + std::to_string(patch.zone) +
                    " face " + bcdef::facedir_to_string(patch.face) +
                    ". Existing: " + existing_desc +
                    " range " + range_to_string(existing.range) +
                    ", New: " + spec->family +
                    " range " + range_to_string(new_range) +
                    ", Overlap: " + range_to_string(overlap));
            } else {
                // Auto BC: skip if any overlap with existing BCs.
                overlap_skip = true;
                break;
            }
        }
        if (overlap_skip)
            continue;

        if (it == specs.end() && (autowall || autofarfield)) {
            bool has_any = false;
            for (const auto& existing : zone_index.entries) {
                if (existing.has_range && ranges_equal(existing.range, new_range)) {
                    has_any = true;
                    break;
                }
            }
            if (has_any)
                continue;
        }

        ensure_family(fn, B, spec->family, spec->type, family_cache);

        int& count = patch_counts[key];
        std::string bc_name;
        bool skip_write = false;
        while (true) {
            ++count;
            bc_name = bc_patch_name(spec->family, patch.face, count);
            auto name_it = zone_index.by_name.find(bc_name);
            if (name_it != zone_index.by_name.end()) {
                const ExistingBC& existing = zone_index.entries[name_it->second];
                if (existing.has_range && ranges_equal(existing.range, new_range)) {
                    skip_write = true;
                    break;
                }
                continue;
            }
            break;
        }
        if (skip_write)
            continue;

        int bc = 0;
        if (mesh.cell_dim() == 2) {
            cgsize_t pnts[4] = {
                new_range[0],
                new_range[1],
                new_range[3],
                new_range[4]
            };
            CG_CALL(cg_boco_write(fn, B, patch.zone, bc_name.c_str(),
                                  CGNS_ENUMV(FamilySpecified),
                                  CGNS_ENUMV(::PointRange),
                                  2, pnts, &bc),
                    "cg_boco_write failed");
        } else {
            cgsize_t pnts[6] = {
                new_range[0],
                new_range[1],
                new_range[2],
                new_range[3],
                new_range[4],
                new_range[5]
            };
            CG_CALL(cg_boco_write(fn, B, patch.zone, bc_name.c_str(),
                                  CGNS_ENUMV(FamilySpecified),
                                  CGNS_ENUMV(::PointRange),
                                  2, pnts, &bc),
                    "cg_boco_write failed");
        }

        CG_CALL(cg_boco_gridlocation_write(fn, B, patch.zone, bc,
                                           face_grid_location(patch.face, mesh.cell_dim())),
                "cg_boco_gridlocation_write failed");

        CG_CALL(cg_goto(fn, B, "Zone_t", patch.zone,
                        "ZoneBC_t", 1, "BC_t", bc, "end"),
                "cg_goto failed (BC_t)");
        CG_CALL(cg_famname_write(spec->family.c_str()),
                "cg_famname_write failed");

        ExistingBC new_entry;
        new_entry.name = bc_name;
        new_entry.family = spec->family;
        new_entry.type = spec->type;
        new_entry.range = new_range;
        new_entry.has_range = true;
        zone_index.by_name[new_entry.name] = zone_index.entries.size();
        zone_index.entries.push_back(new_entry);
    }
}

} // namespace bcdef::boundary
