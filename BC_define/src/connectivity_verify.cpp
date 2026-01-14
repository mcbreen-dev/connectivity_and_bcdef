//─────────────────────────────────────────────────────────────
// File: src/connectivity_verify.cpp
// Match verification and transform resolution
//─────────────────────────────────────────────────────────────
#include "connectivity_internal.hpp"
#include "logger.hpp"

#include <algorithm>
#include <limits>
#include <thread>

namespace fs::conn {

bool derive_uv_mapping(const std::vector<UVMatch>& matches,
                       bool need_u,
                       bool need_v,
                       int& mu_u,
                       int& mu_v,
                       int& mv_u,
                       int& mv_v)
{
    if (matches.empty())
        return false;

    const UVMatch& m0 = matches.front();
    mu_u = 0; mu_v = 0; mv_u = 0; mv_v = 0;

    auto sgn = [](int v) { return v >= 0 ? 1 : -1; };

    for (size_t i = 1; i < matches.size(); ++i) {
        const UVMatch& m = matches[i];
        int du = m.au - m0.au;
        int dv = m.av - m0.av;
        int dbu = m.bu - m0.bu;
        int dbv = m.bv - m0.bv;

        if (need_u && (mu_u == 0 && mu_v == 0)) {
            if (dv == 0 && du != 0) {
                if (dbv == 0 && dbu != 0) {
                    mu_u = sgn(dbu) * sgn(du);
                    mu_v = 0;
                } else if (dbu == 0 && dbv != 0) {
                    mu_u = 0;
                    mu_v = sgn(dbv) * sgn(du);
                }
            }
        }

        if (need_v && (mv_u == 0 && mv_v == 0)) {
            if (du == 0 && dv != 0) {
                if (dbv == 0 && dbu != 0) {
                    mv_u = sgn(dbu) * sgn(dv);
                    mv_v = 0;
                } else if (dbu == 0 && dbv != 0) {
                    mv_u = 0;
                    mv_v = sgn(dbv) * sgn(dv);
                }
            }
        }

        if ((!need_u || (mu_u != 0 || mu_v != 0)) &&
            (!need_v || (mv_u != 0 || mv_v != 0)))
            break;
    }

    if (need_u && (mu_u == 0 && mu_v == 0))
        return false;
    if (need_v && (mv_u == 0 && mv_v == 0))
        return false;

    if ((mu_u != 0 || mu_v != 0) && (mv_u != 0 || mv_v != 0)) {
        if (mu_u == mv_u && mu_v == mv_v)
            return false;
        if (mu_u * mv_u + mu_v * mv_v != 0)
            return false;
    }

    return true;
}

int find_match_in_face(const FaceHash& hash,
                       const FaceData& faces,
                       double tx,
                       double ty,
                       double tz,
                       double tol2,
                       double cell)
{
    CellKey key = cell_key(tx, ty, tz, cell);
    double best = std::numeric_limits<double>::max();
    int best_idx = -1;

    for (int dx = -1; dx <= 1; ++dx) {
        for (int dy = -1; dy <= 1; ++dy) {
            for (int dz = -1; dz <= 1; ++dz) {
                CellKey probe{key.ix + dx, key.iy + dy, key.iz + dz};
                auto it = hash.cells.find(probe);
                if (it == hash.cells.end())
                    continue;
                for (int candidate : it->second) {
                    double d = distance_sq(faces.x[candidate], faces.y[candidate], faces.z[candidate],
                                           tx, ty, tz);
                    if (d < best || (d == best && candidate < best_idx)) {
                        best = d;
                        best_idx = candidate;
                    }
                }
            }
        }
    }

    if (best_idx != -1 && best <= tol2)
        return best_idx;
    return -1;
}

std::vector<MatchPair> collect_match_pairs(const FaceBuildResult& build,
                                           const std::vector<FaceSummary>& summaries,
                                           const HashBuildResult& hash,
                                           double tol,
                                           int threads)
{
    std::vector<MatchPair> all_matches;
    if (build.faces.empty())
        return all_matches;

    double tol2 = tol * tol;
    double hash_cell = hash.hash_cell;
    const auto& pairs = hash.pairs;

    auto collect_matches_for_pair = [&](int faceA, int faceB,
                                        std::vector<MatchPair>& out) {
        size_t startA = build.cell_start[static_cast<size_t>(faceA)];
        int countA = build.cell_count[static_cast<size_t>(faceA)];
        size_t startB = build.cell_start[static_cast<size_t>(faceB)];
        int countB = build.cell_count[static_cast<size_t>(faceB)];
        if (countA == 0 || countB == 0)
            return;

        const FaceSummary& A = summaries[static_cast<size_t>(faceA)];
        const FaceSummary& B = summaries[static_cast<size_t>(faceB)];
        const FaceHash& hashB = hash.face_hashes[static_cast<size_t>(faceB)];

        std::vector<UVMatch> matches;
        matches.reserve(8);
        bool need_u = A.uLen > 1;
        bool need_v = A.vLen > 1;
        int mu_u = 0, mu_v = 0, mv_u = 0, mv_v = 0;
        bool mapping_ok = false;

        for (int ii = 0; ii < countA; ++ii) {
            int idx = build.face_cells[startA + static_cast<size_t>(ii)];
            int b_idx = find_match_in_face(hashB, build.faces,
                                           build.faces.x[idx], build.faces.y[idx], build.faces.z[idx],
                                           tol2, hash_cell);
            if (b_idx < 0)
                continue;
            matches.push_back(UVMatch{
                build.faces.u[idx], build.faces.v[idx],
                build.faces.u[b_idx], build.faces.v[b_idx], idx, b_idx});
            if (!mapping_ok && matches.size() >= 3) {
                mapping_ok = derive_uv_mapping(matches, need_u, need_v, mu_u, mu_v, mv_u, mv_v);
                if (mapping_ok)
                    break;
            }
        }

        if (mapping_ok) {
            const UVMatch& origin = matches.front();
            const FaceGrid& gridB = build.grids[static_cast<size_t>(faceB)];
            std::vector<MatchPair> local;
            local.reserve(static_cast<size_t>(countA));
            for (int ii = 0; ii < countA; ++ii) {
                int idx = build.face_cells[startA + static_cast<size_t>(ii)];
                int du = build.faces.u[idx] - origin.au;
                int dv = build.faces.v[idx] - origin.av;
                int bu = origin.bu + du * mu_u + dv * mv_u;
                int bv = origin.bv + du * mu_v + dv * mv_v;
                if (bu < 0 || bu >= B.uLen || bv < 0 || bv >= B.vLen)
                    continue;
                int b_idx = gridB.index[face_grid_index(gridB, bu, bv)];
                if (b_idx < 0)
                    continue;
                if (distance_sq(build.faces.x[idx], build.faces.y[idx], build.faces.z[idx],
                                build.faces.x[b_idx], build.faces.y[b_idx], build.faces.z[b_idx]) <= tol2)
                    local.push_back(MatchPair{idx, b_idx});
            }
            if (!local.empty()) {
                out.insert(out.end(), local.begin(), local.end());
                return;
            }
        }

        for (int ii = 0; ii < countA; ++ii) {
            int idx = build.face_cells[startA + static_cast<size_t>(ii)];
            int b_idx = find_match_in_face(hashB, build.faces,
                                           build.faces.x[idx], build.faces.y[idx], build.faces.z[idx],
                                           tol2, hash_cell);
            if (b_idx >= 0)
                out.push_back(MatchPair{idx, b_idx});
        }
    };

    auto process_pairs = [&](size_t begin, size_t end, std::vector<MatchPair>& out) {
        for (size_t i = begin; i < end; ++i) {
            int faceA = pairs[i].a;
            int faceB = pairs[i].b;
            bool swapped = false;
            if (build.cell_count[static_cast<size_t>(faceA)] >
                build.cell_count[static_cast<size_t>(faceB)]) {
                std::swap(faceA, faceB);
                swapped = true;
            }
            std::vector<MatchPair> local;
            collect_matches_for_pair(faceA, faceB, local);
            if (swapped) {
                for (const auto& m : local)
                    out.push_back(MatchPair{m.b, m.a});
            } else {
                out.insert(out.end(), local.begin(), local.end());
            }
        }
    };

    int thread_count = threads;
    if (thread_count <= 0)
        thread_count = static_cast<int>(std::thread::hardware_concurrency());
    if (thread_count < 1)
        thread_count = 1;

    std::vector<std::vector<MatchPair>> thread_matches(static_cast<size_t>(thread_count));
    size_t pair_total = pairs.size();
    if (thread_count == 1 || pair_total == 0) {
        process_pairs(0, pair_total, thread_matches[0]);
    } else {
        std::vector<std::thread> workers;
        workers.reserve(static_cast<size_t>(thread_count));
        size_t chunk = (pair_total + static_cast<size_t>(thread_count) - 1) /
                       static_cast<size_t>(thread_count);
        for (int t = 0; t < thread_count; ++t) {
            size_t begin = static_cast<size_t>(t) * chunk;
            if (begin >= pair_total)
                break;
            size_t end = std::min(pair_total, begin + chunk);
            workers.emplace_back(process_pairs, begin, end,
                                 std::ref(thread_matches[static_cast<size_t>(t)]));
        }
        for (auto& worker : workers)
            worker.join();
    }

    size_t total_matches = 0;
    for (const auto& list : thread_matches)
        total_matches += list.size();
    all_matches.reserve(total_matches);
    for (auto& list : thread_matches)
        all_matches.insert(all_matches.end(), list.begin(), list.end());

    std::sort(all_matches.begin(), all_matches.end(),
              [](const MatchPair& lhs, const MatchPair& rhs) {
                  if (lhs.a != rhs.a) return lhs.a < rhs.a;
                  return lhs.b < rhs.b;
              });

    return all_matches;
}

void resolve_matches_cgns(const Mesh& mesh, Logger& log,
                          FaceBuildResult& build,
                          const std::vector<MatchPair>& matches,
                          std::vector<BoundaryRecord>& boundary_records,
                          std::vector<ConnPatch>& conn_patches)
{
    std::vector<NeighborInfo> otherInfo(build.faces.size());

    for (const auto& match : matches) {
        if (otherInfo[static_cast<size_t>(match.a)].block != -1 ||
            otherInfo[static_cast<size_t>(match.b)].block != -1)
            continue;
        otherInfo[match.a] = NeighborInfo{
            build.faces.block[match.b], build.faces.face[match.b],
            build.faces.i[match.b], build.faces.j[match.b], build.faces.k[match.b]};
        otherInfo[match.b] = NeighborInfo{
            build.faces.block[match.a], build.faces.face[match.a],
            build.faces.i[match.a], build.faces.j[match.a], build.faces.k[match.a]};
    }

    const int nBlocks = static_cast<int>(mesh.zones().size());
    for (int b = 0; b < nBlocks; ++b)
    {
        const Zone& z = mesh.zones()[b];
        for (int face = 0; face < 6; ++face)
        {
            int face_id = b * 6 + face;
            FaceGrid& grid = build.grids[static_cast<size_t>(face_id)];
            if (grid.index.empty()) continue;

            int uLen = grid.uLen;
            int vLen = grid.vLen;

            int masterCount = 0;
            int uStart = 0;
            int vStart = 0;
            int uEnd = 0;
            int vEnd = 0;
            bool complete = false;

            while (!complete) {
                const auto startIdx = grid.index[face_grid_index(grid, uStart, vStart)];
                NeighborInfo curOther = otherInfo[startIdx];
                NeighborInfo startOther = curOther;

                int iDirIndex = 0;
                int jDirIndex = 0;

                while (uEnd < uLen - 1) {
                    int nextU = uEnd + 1;
                    const auto idx = grid.index[face_grid_index(grid, nextU, vStart)];
                    NeighborInfo nextOther = otherInfo[idx];

                    int delI = std::abs(nextOther.i - curOther.i);
                    int delJ = std::abs(nextOther.j - curOther.j);
                    int delK = std::abs(nextOther.k - curOther.k);

                    bool ok = false;
                    if (nextOther.block == curOther.block && nextOther.block == -1) {
                        ok = true;
                    } else if (nextOther.block == curOther.block &&
                               nextOther.face == curOther.face &&
                               (delI == 1 || delJ == 1 || delK == 1)) {
                        ok = true;
                    }

                    if (!ok) break;

                    if (iDirIndex == 0 && nextOther.block != -1) {
                        if (delI == 1) iDirIndex = signed_axis(1, nextOther.i - curOther.i);
                        else if (delJ == 1) iDirIndex = signed_axis(2, nextOther.j - curOther.j);
                        else if (delK == 1) iDirIndex = signed_axis(3, nextOther.k - curOther.k);
                        else throw std::runtime_error("Unexpected I-direction delta");
                    }

                    curOther = nextOther;
                    uEnd = nextU;
                }

                curOther = otherInfo[grid.index[face_grid_index(grid, uEnd, vStart)]];

                while (vEnd < vLen - 1) {
                    int nextV = vEnd + 1;
                    const auto idx = grid.index[face_grid_index(grid, uEnd, nextV)];
                    NeighborInfo nextOther = otherInfo[idx];

                    int delI = std::abs(nextOther.i - curOther.i);
                    int delJ = std::abs(nextOther.j - curOther.j);
                    int delK = std::abs(nextOther.k - curOther.k);

                    bool ok = false;
                    if (nextOther.block == curOther.block && nextOther.block == -1) {
                        ok = true;
                    } else if (nextOther.block == curOther.block &&
                               nextOther.face == curOther.face &&
                               (delI == 1 || delJ == 1 || delK == 1)) {
                        ok = true;
                    }

                    if (!ok) break;

                    if (jDirIndex == 0 && nextOther.block != -1) {
                        if (delI == 1) jDirIndex = signed_axis(1, nextOther.i - curOther.i);
                        else if (delJ == 1) jDirIndex = signed_axis(2, nextOther.j - curOther.j);
                        else if (delK == 1) jDirIndex = signed_axis(3, nextOther.k - curOther.k);
                        else throw std::runtime_error("Unexpected J-direction delta");
                    }

                    curOther = nextOther;
                    vEnd = nextV;
                }

                for (int uu = uStart; uu <= uEnd; ++uu)
                    for (int vv = vStart; vv <= vEnd; ++vv)
                        grid.consumed[face_grid_index(grid, uu, vv)] = 1;

                IJK recvBegin{};
                IJK recvEnd{};
                fill_point_range(z, face, uStart, vStart, uEnd, vEnd, recvBegin, recvEnd);

                if (z.nk() == 1) {
                    recvBegin[2] = 1;
                    recvEnd[2] = 1;
                }

                if (curOther.block == -1) {
                    bool isFull = (uStart == 0 && vStart == 0 && uEnd == uLen - 1 && vEnd == vLen - 1);
                    BoundaryRecord entry;
                    entry.zone = b + 1;
                    entry.face = face_dir_from_id(face);
                    entry.full = isFull;
                    entry.begin = recvBegin;
                    entry.end = recvEnd;
                    boundary_records.push_back(entry);
                } else {
                    int donorFace = curOther.face;

                    if (iDirIndex == 0 && jDirIndex != 0) {
                        if (is_i_face(donorFace)) {
                            iDirIndex = (std::abs(jDirIndex) == 2) ? 3 : 2;
                        } else if (is_j_face(donorFace)) {
                            iDirIndex = (std::abs(jDirIndex) == 1) ? 3 : 1;
                        } else {
                            iDirIndex = (std::abs(jDirIndex) == 1) ? 2 : 1;
                        }
                    } else if (iDirIndex != 0 && jDirIndex == 0) {
                        if (is_i_face(donorFace)) {
                            jDirIndex = (std::abs(iDirIndex) == 2) ? 3 : 2;
                        } else if (is_j_face(donorFace)) {
                            jDirIndex = (std::abs(iDirIndex) == 1) ? 3 : 1;
                        } else {
                            jDirIndex = (std::abs(iDirIndex) == 1) ? 2 : 1;
                        }
                    } else if (iDirIndex == 0 || jDirIndex == 0) {
                        if (is_i_face(donorFace)) {
                            iDirIndex = 2;
                            jDirIndex = 3;
                        } else if (is_j_face(donorFace)) {
                            iDirIndex = 1;
                            jDirIndex = 3;
                        } else {
                            iDirIndex = 1;
                            jDirIndex = 2;
                        }
                    }

                    int transform[3] = {0, 0, 0};
                    if (is_i_face(face)) {
                        transform[1] = iDirIndex;
                        transform[2] = jDirIndex;
                    } else if (is_j_face(face)) {
                        transform[0] = iDirIndex;
                        transform[2] = jDirIndex;
                    } else {
                        transform[0] = iDirIndex;
                        transform[1] = jDirIndex;
                    }
                    set_normal_transform(transform, face, donorFace);

                    IJK donorBegin{};
                    IJK donorEnd{};
                    IJK donorStart = {startOther.i, startOther.j, startOther.k};
                    compute_donor_range(transform, recvBegin, recvEnd, donorStart, donorBegin, donorEnd);

                    const Zone& donorZ = mesh.zones()[curOther.block];
                    if (donorZ.nk() == 1) {
                        donorBegin[2] = 1;
                        donorEnd[2] = 1;
                    }

                    int recvZone = b + 1;
                    int donorZone = curOther.block + 1;
                    if (should_write(recvZone, donorZone, face, donorFace)) {
                        int transform_cgns[3] = {transform[0], transform[1], transform[2]};
                        if (z.nk() == 1 || donorZ.nk() == 1) {
                            transform_cgns[2] = 0;
                        }
                        bool isFull = (uStart == 0 && vStart == 0 && uEnd == uLen - 1 && vEnd == vLen - 1);
                        std::string recvFaceStr = facedir_to_string(face_dir_from_id(face));
                        std::string donorFaceStr = facedir_to_string(face_dir_from_id(donorFace));

                        if (isFull) {
                            log.info("interior Z" + std::to_string(recvZone) + " " + recvFaceStr +
                                     " ↔ Z" + std::to_string(donorZone) + " " + donorFaceStr +
                                     "  (T=" + std::to_string(transform_cgns[0]) + "," +
                                     std::to_string(transform_cgns[1]) + "," +
                                     std::to_string(transform_cgns[2]) + ")");
                        } else {
                            log.info("interior PATCH Z" + std::to_string(recvZone) + " " + recvFaceStr +
                                     " ↔ Z" + std::to_string(donorZone) + " " + donorFaceStr +
                                     "  (T=" + std::to_string(transform_cgns[0]) + "," +
                                     std::to_string(transform_cgns[1]) + "," +
                                     std::to_string(transform_cgns[2]) + ")");
                        }

                        ConnPatch cp;
                        cp.recvZone = recvZone;
                        cp.donorZone = donorZone;
                        cp.recvDir = face_dir_from_id(face);
                        cp.donorDir = face_dir_from_id(donorFace);
                        cp.offU = uStart;
                        cp.offV = vStart;
                        cp.sizeU = (uEnd - uStart + 2);
                        cp.sizeV = (vEnd - vStart + 2);
                        cp.ori = 0;
                        cp.transform[0] = transform_cgns[0];
                        cp.transform[1] = transform_cgns[1];
                        cp.transform[2] = transform_cgns[2];
                        cp.recvRange.begin = recvBegin;
                        cp.recvRange.end = recvEnd;
                        cp.donorRange.begin = donorBegin;
                        cp.donorRange.end = donorEnd;
                        cp.useExplicitRanges = true;

                        conn_patches.push_back(cp);
                    }
                }

                while (masterCount < uLen * vLen) {
                    masterCount += 1;
                    int uu = (masterCount - 1) % uLen;
                    int vv = (masterCount - 1) / uLen;

                    if (!grid.consumed[face_grid_index(grid, uu, vv)]) {
                        uStart = uu;
                        vStart = vv;
                        uEnd = uStart;
                        vEnd = vStart;
                        break;
                    }
                }

                if (masterCount == uLen * vLen) {
                    complete = true;
                }
            }
        }
    }
}

void resolve_matches_plot3d(const std::vector<Plot3DZone>& zones, Logger& log,
                            FaceBuildResult& build,
                            const std::vector<MatchPair>& matches,
                            std::vector<BoundaryRecord>& boundary_records,
                            std::vector<ConnPatch>& conn_patches)
{
    std::vector<NeighborInfo> otherInfo(build.faces.size());

    for (const auto& match : matches) {
        if (otherInfo[static_cast<size_t>(match.a)].block != -1 ||
            otherInfo[static_cast<size_t>(match.b)].block != -1)
            continue;
        otherInfo[match.a] = NeighborInfo{
            build.faces.block[match.b], build.faces.face[match.b],
            build.faces.i[match.b], build.faces.j[match.b], build.faces.k[match.b]};
        otherInfo[match.b] = NeighborInfo{
            build.faces.block[match.a], build.faces.face[match.a],
            build.faces.i[match.a], build.faces.j[match.a], build.faces.k[match.a]};
    }

    const int nBlocks = static_cast<int>(zones.size());
    for (int b = 0; b < nBlocks; ++b)
    {
        const Plot3DZone& z = zones[b];
        for (int face = 0; face < 6; ++face)
        {
            int face_id = b * 6 + face;
            FaceGrid& grid = build.grids[static_cast<size_t>(face_id)];
            if (grid.index.empty()) continue;

            int uLen = grid.uLen;
            int vLen = grid.vLen;

            int masterCount = 0;
            int uStart = 0;
            int vStart = 0;
            int uEnd = 0;
            int vEnd = 0;
            bool complete = false;

            while (!complete) {
                const auto startIdx = grid.index[face_grid_index(grid, uStart, vStart)];
                NeighborInfo curOther = otherInfo[startIdx];
                NeighborInfo startOther = curOther;

                int iDirIndex = 0;
                int jDirIndex = 0;

                while (uEnd < uLen - 1) {
                    int nextU = uEnd + 1;
                    const auto idx = grid.index[face_grid_index(grid, nextU, vStart)];
                    NeighborInfo nextOther = otherInfo[idx];

                    int delI = std::abs(nextOther.i - curOther.i);
                    int delJ = std::abs(nextOther.j - curOther.j);
                    int delK = std::abs(nextOther.k - curOther.k);

                    bool ok = false;
                    if (nextOther.block == curOther.block && nextOther.block == -1) {
                        ok = true;
                    } else if (nextOther.block == curOther.block &&
                               nextOther.face == curOther.face &&
                               (delI == 1 || delJ == 1 || delK == 1)) {
                        ok = true;
                    }

                    if (!ok) break;

                    if (iDirIndex == 0 && nextOther.block != -1) {
                        if (delI == 1) iDirIndex = signed_axis(1, nextOther.i - curOther.i);
                        else if (delJ == 1) iDirIndex = signed_axis(2, nextOther.j - curOther.j);
                        else if (delK == 1) iDirIndex = signed_axis(3, nextOther.k - curOther.k);
                        else throw std::runtime_error("Unexpected I-direction delta");
                    }

                    curOther = nextOther;
                    uEnd = nextU;
                }

                curOther = otherInfo[grid.index[face_grid_index(grid, uEnd, vStart)]];

                while (vEnd < vLen - 1) {
                    int nextV = vEnd + 1;
                    const auto idx = grid.index[face_grid_index(grid, uEnd, nextV)];
                    NeighborInfo nextOther = otherInfo[idx];

                    int delI = std::abs(nextOther.i - curOther.i);
                    int delJ = std::abs(nextOther.j - curOther.j);
                    int delK = std::abs(nextOther.k - curOther.k);

                    bool ok = false;
                    if (nextOther.block == curOther.block && nextOther.block == -1) {
                        ok = true;
                    } else if (nextOther.block == curOther.block &&
                               nextOther.face == curOther.face &&
                               (delI == 1 || delJ == 1 || delK == 1)) {
                        ok = true;
                    }

                    if (!ok) break;

                    if (jDirIndex == 0 && nextOther.block != -1) {
                        if (delI == 1) jDirIndex = signed_axis(1, nextOther.i - curOther.i);
                        else if (delJ == 1) jDirIndex = signed_axis(2, nextOther.j - curOther.j);
                        else if (delK == 1) jDirIndex = signed_axis(3, nextOther.k - curOther.k);
                        else throw std::runtime_error("Unexpected J-direction delta");
                    }

                    curOther = nextOther;
                    vEnd = nextV;
                }

                for (int uu = uStart; uu <= uEnd; ++uu)
                    for (int vv = vStart; vv <= vEnd; ++vv)
                        grid.consumed[face_grid_index(grid, uu, vv)] = 1;

                IJK recvBegin{};
                IJK recvEnd{};
                fill_point_range(z, face, uStart, vStart, uEnd, vEnd, recvBegin, recvEnd);

                if (z.nk() == 1) {
                    recvBegin[2] = 1;
                    recvEnd[2] = 1;
                }

                if (curOther.block == -1) {
                    bool isFull = (uStart == 0 && vStart == 0 && uEnd == uLen - 1 && vEnd == vLen - 1);
                    BoundaryRecord entry;
                    entry.zone = b + 1;
                    entry.face = face_dir_from_id(face);
                    entry.full = isFull;
                    entry.begin = recvBegin;
                    entry.end = recvEnd;
                    boundary_records.push_back(entry);
                } else {
                    int donorFace = curOther.face;

                    if (iDirIndex == 0 && jDirIndex != 0) {
                        if (is_i_face(donorFace)) {
                            iDirIndex = (std::abs(jDirIndex) == 2) ? 3 : 2;
                        } else if (is_j_face(donorFace)) {
                            iDirIndex = (std::abs(jDirIndex) == 1) ? 3 : 1;
                        } else {
                            iDirIndex = (std::abs(jDirIndex) == 1) ? 2 : 1;
                        }
                    } else if (iDirIndex != 0 && jDirIndex == 0) {
                        if (is_i_face(donorFace)) {
                            jDirIndex = (std::abs(iDirIndex) == 2) ? 3 : 2;
                        } else if (is_j_face(donorFace)) {
                            jDirIndex = (std::abs(iDirIndex) == 1) ? 3 : 1;
                        } else {
                            jDirIndex = (std::abs(iDirIndex) == 1) ? 2 : 1;
                        }
                    } else if (iDirIndex == 0 || jDirIndex == 0) {
                        if (is_i_face(donorFace)) {
                            iDirIndex = 2;
                            jDirIndex = 3;
                        } else if (is_j_face(donorFace)) {
                            iDirIndex = 1;
                            jDirIndex = 3;
                        } else {
                            iDirIndex = 1;
                            jDirIndex = 2;
                        }
                    }

                    int transform[3] = {0, 0, 0};
                    if (is_i_face(face)) {
                        transform[1] = iDirIndex;
                        transform[2] = jDirIndex;
                    } else if (is_j_face(face)) {
                        transform[0] = iDirIndex;
                        transform[2] = jDirIndex;
                    } else {
                        transform[0] = iDirIndex;
                        transform[1] = jDirIndex;
                    }
                    set_normal_transform(transform, face, donorFace);

                    IJK donorBegin{};
                    IJK donorEnd{};
                    IJK donorStart = {startOther.i, startOther.j, startOther.k};
                    compute_donor_range(transform, recvBegin, recvEnd, donorStart, donorBegin, donorEnd);

                    const Plot3DZone& donorZ = zones[static_cast<size_t>(curOther.block)];
                    if (donorZ.nk() == 1) {
                        donorBegin[2] = 1;
                        donorEnd[2] = 1;
                    }

                    int recvZone = b + 1;
                    int donorZone = curOther.block + 1;
                    if (should_write(recvZone, donorZone, face, donorFace)) {
                        int transform_cgns[3] = {transform[0], transform[1], transform[2]};
                        if (z.nk() == 1 || donorZ.nk() == 1) {
                            transform_cgns[2] = 0;
                        }
                        bool isFull = (uStart == 0 && vStart == 0 && uEnd == uLen - 1 && vEnd == vLen - 1);
                        std::string recvFaceStr = facedir_to_string(face_dir_from_id(face));
                        std::string donorFaceStr = facedir_to_string(face_dir_from_id(donorFace));

                        if (isFull) {
                            log.info("interior Z" + std::to_string(recvZone) + " " + recvFaceStr +
                                     " ↔ Z" + std::to_string(donorZone) + " " + donorFaceStr +
                                     "  (T=" + std::to_string(transform_cgns[0]) + "," +
                                     std::to_string(transform_cgns[1]) + "," +
                                     std::to_string(transform_cgns[2]) + ")");
                        } else {
                            log.info("interior PATCH Z" + std::to_string(recvZone) + " " + recvFaceStr +
                                     " ↔ Z" + std::to_string(donorZone) + " " + donorFaceStr +
                                     "  (T=" + std::to_string(transform_cgns[0]) + "," +
                                     std::to_string(transform_cgns[1]) + "," +
                                     std::to_string(transform_cgns[2]) + ")");
                        }

                        ConnPatch cp;
                        cp.recvZone = recvZone;
                        cp.donorZone = donorZone;
                        cp.recvDir = face_dir_from_id(face);
                        cp.donorDir = face_dir_from_id(donorFace);
                        cp.offU = uStart;
                        cp.offV = vStart;
                        cp.sizeU = (uEnd - uStart + 2);
                        cp.sizeV = (vEnd - vStart + 2);
                        cp.ori = 0;
                        cp.transform[0] = transform_cgns[0];
                        cp.transform[1] = transform_cgns[1];
                        cp.transform[2] = transform_cgns[2];
                        cp.recvRange.begin = recvBegin;
                        cp.recvRange.end = recvEnd;
                        cp.donorRange.begin = donorBegin;
                        cp.donorRange.end = donorEnd;
                        cp.useExplicitRanges = true;

                        conn_patches.push_back(cp);
                    }
                }

                while (masterCount < uLen * vLen) {
                    masterCount += 1;
                    int uu = (masterCount - 1) % uLen;
                    int vv = (masterCount - 1) / uLen;

                    if (!grid.consumed[face_grid_index(grid, uu, vv)]) {
                        uStart = uu;
                        vStart = vv;
                        uEnd = uStart;
                        vEnd = vStart;
                        break;
                    }
                }

                if (masterCount == uLen * vLen) {
                    complete = true;
                }
            }
        }
    }
}

} // namespace fs::conn
