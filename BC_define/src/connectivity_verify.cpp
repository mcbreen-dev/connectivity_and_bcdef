//─────────────────────────────────────────────────────────────
// File: src/connectivity_verify.cpp
// Match verification and transform resolution
//─────────────────────────────────────────────────────────────
//
// This file performs two jobs that sit after candidate generation:
//
// (1) Collecting detailed point-to-point matches between two candidate
//     faces using a spatial hash (FaceHash) built in connectivity_hash.cpp.
//
// (2) Turning those per-point matches into:
//       - BoundaryRecord entries for unmatched regions (BC targets)
//       - ConnPatch entries for matched regions (1-to-1 connectivity targets)
//
// Data flow for CGNS:
//
//   build_faces_cgns()                    (connectivity_face_build.cpp)
//       -> FaceBuildResult build
//   build_face_summaries()                (connectivity_face_build.cpp)
//       -> per-face bounds, centroids, normals
//   build_hash_candidates()               (connectivity_hash.cpp)
//       -> per-face point hash + candidate face pairs
//   collect_match_pairs()                 (this file)
//       -> flat list of matched FaceData index pairs
//   resolve_matches_cgns()                (this file)
//       -> merges match information into:
//            * BoundaryRecord (unmatched face patches)
//            * ConnPatch (interior connectivity patches)
//
// Plot3D path is identical except that face extraction and handedness
// checks use Plot3DZone and build_faces_plot3d().
//
// Indexing and sizes:
//
//   - Face IDs are "block*6 + face", where face is 0..5.
//   - FaceGrid is sized uLen*vLen with 0-based (u,v) indexing.
//   - FaceData is a struct-of-arrays (SoA) inside FaceBuildResult::faces.
//     Each FaceData entry represents one face cell center point.
//   - build.face_cells is a flat list of FaceData indices grouped by face,
//     with per-face start/count stored in build.cell_start/cell_count.
//   - NeighborInfo stores a *single* neighbor per FaceData index and is
//     initialized to block=-1 meaning “no neighbor” (boundary).
//
// Patch semantics:
//
//   - Resolution groups a face into rectangular patches in (u,v) index space.
//   - A patch is labeled as:
//       * boundary if its points have no neighbor (block == -1)
//       * interior if its points all connect to one donor block+face with
//         consistent step directions.
//
// 2D zones (nk==1):
//
//   - Only faces 0..3 are populated by face_dims() / build_faces_*.
//   - vLen is 1 on those populated faces.
//   - transform[2] is forced to 0 when either side is 2D, and donor/recv
//     k indices are pinned to 1.
//
//─────────────────────────────────────────────────────────────
#include "connectivity_internal.hpp"
#include "logger.hpp"

#include <algorithm>
#include <limits>
#include <thread>

namespace fs::conn {

/*=====================================================================
// derive_uv_mapping
//
// Attempt to infer a linear mapping between face-local integer grid
// coordinates:
//
//   Given matches between points on Face A and Face B:
//
//       (au,av)  -> (bu,bv)
//
// Infer how A's u and v step directions map onto B's u and v axes:
//
//   bu = b0u + (au-a0u) * mu_u + (av-a0v) * mv_u
//   bv = b0v + (au-a0u) * mu_v + (av-a0v) * mv_v
//
// Where:
//   - (mu_u, mu_v) is the image of +u_A in B coordinates.
//   - (mv_u, mv_v) is the image of +v_A in B coordinates.
//   - Each vector is axis-aligned (one component is 0, the other is ±1),
//     and mu and mv must be orthogonal (dot == 0) and not identical.
//
// Parameters:
//   matches : small set of point matches (UVMatch) for one face pair
//   need_u  : whether u direction needs resolution (A.uLen > 1)
//   need_v  : whether v direction needs resolution (A.vLen > 1)
//   outputs : mapping coefficients (mu_u, mu_v, mv_u, mv_v)
//
// Return value:
//   true  -> mapping inferred and valid
//   false -> mapping could not be inferred or is inconsistent
//
// Notes on minimum evidence:
//   - The routine requires observing at least one match where du!=0,dv==0
//     (pure u step in A) to infer mu_*.
//   - Similarly, requires one match where du==0,dv!=0 to infer mv_*.
//   - In 2D faces (vLen==1), need_v is false and only u direction is required.
//=====================================================================*/
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

        // Differences in face-local coordinates on A and B relative to origin m0.
        int du = m.au - m0.au;
        int dv = m.av - m0.av;
        int dbu = m.bu - m0.bu;
        int dbv = m.bv - m0.bv;

        /*-------------------------------------------------------------
          Infer mapping for A's u direction when a “pure u” step on A
          is observed (dv==0, du!=0).

          If B also stepped purely in bu or bv, record the mapping.
        -------------------------------------------------------------*/
        if (need_u && (mu_u == 0 && mu_v == 0)) {
            if (dv == 0 && du != 0) {
                if (dbv == 0 && dbu != 0) {
                    // +u_A corresponds to ±u_B
                    mu_u = sgn(dbu) * sgn(du);
                    mu_v = 0;
                } else if (dbu == 0 && dbv != 0) {
                    // +u_A corresponds to ±v_B
                    mu_u = 0;
                    mu_v = sgn(dbv) * sgn(du);
                }
            }
        }

        /*-------------------------------------------------------------
          Infer mapping for A's v direction when a “pure v” step on A
          is observed (du==0, dv!=0).
        -------------------------------------------------------------*/
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

        // Stop once all required directions are resolved.
        if ((!need_u || (mu_u != 0 || mu_v != 0)) &&
            (!need_v || (mv_u != 0 || mv_v != 0)))
            break;
    }

    // Required direction(s) could not be inferred.
    if (need_u && (mu_u == 0 && mu_v == 0))
        return false;
    if (need_v && (mv_u == 0 && mv_v == 0))
        return false;

    /*-------------------------------------------------------------
      If both directions exist, enforce a valid 2D basis in B:
        - mu and mv must not be identical
        - mu · mv must be 0 (orthogonal)
    -------------------------------------------------------------*/
    if ((mu_u != 0 || mu_v != 0) && (mv_u != 0 || mv_v != 0)) {
        if (mu_u == mv_u && mu_v == mv_v)
            return false;
        if (mu_u * mv_u + mu_v * mv_v != 0)
            return false;
    }

    return true;
}

/*=====================================================================
// find_match_in_face
//
// Find the nearest point in Face B to a query point (tx,ty,tz) using:
//
//   - A per-face spatial hash: FaceHash::cells
//   - A 3x3x3 neighborhood of hash cells around the query cell
//
// Input:
//   hash  : FaceHash for candidate face B (CellKey -> vector<int>)
//   faces : FaceData SoA holding coordinates for all face points
//   tx,ty,tz : query point coordinates (double)
//   tol2  : squared tolerance (tol*tol)
//   cell  : hash cell size (same used to build FaceHash)
//
// Output:
//   returns FaceData index in B that is closest to query within tol2,
//   or -1 if no candidate lies within tolerance.
//
// Tie-breaking:
//   - If two candidates have equal squared distance, the smaller FaceData
//     index is selected (stable selection).
//=====================================================================*/
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

    // Probe the query cell and all immediate neighbors.
    for (int dx = -1; dx <= 1; ++dx) {
        for (int dy = -1; dy <= 1; ++dy) {
            for (int dz = -1; dz <= 1; ++dz) {
                CellKey probe{key.ix + dx, key.iy + dy, key.iz + dz};
                auto it = hash.cells.find(probe);
                if (it == hash.cells.end())
                    continue;

                // Each bucket stores FaceData indices that quantize into that cell.
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

/*=====================================================================
// collect_match_pairs
//
// For each candidate face pair (A,B) from HashBuildResult::pairs, find
// point-to-point matches between their face-cell centers.
//
// Outputs:
//   - A flat vector<MatchPair> where each entry is (FaceData index a, FaceData index b).
//
// Matching strategy per face pair:
//
//   1) Attempt to infer a (u,v) mapping between the two faces using a small
//      number of nearest-neighbor matches:
//
//         derive_uv_mapping(matches, need_u, need_v, ...)
//
//      If a mapping is found, use it to predict B(u,v) indices for every
//      point on A using FaceGrid::index (O(N)) rather than doing a hash lookup
//      for each point.
//
//      The predicted match is accepted if the squared distance is <= tol^2.
//
//   2) If no mapping is found (or mapping yields no matches), fall back to
//      per-point spatial hash nearest-neighbor search on B.
//
// Work partitioning:
//
//   - The function splits candidate face pairs across threads.
//   - For a given candidate pair (a,b), the smaller face (fewer cells) is chosen
//     as the “query” side to reduce lookup work; MatchPairs are swapped back
//     to preserve a->b meaning.
//
// Ordering:
//
//   - The final output is sorted by (a,b). No explicit dedup is required here
//     because each A point produces at most one match in this routine.
//=====================================================================*/
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

    /*-------------------------------------------------------------
      Collect matches for a single face pair (faceA, faceB).

      Output is appended to 'out' (local per-thread vector) as MatchPair{a_idx, b_idx}.
    -------------------------------------------------------------*/
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

        // Collect a small set of UVMatch entries as evidence for mapping inference.
        std::vector<UVMatch> matches;
        matches.reserve(8);
        bool need_u = A.uLen > 1;
        bool need_v = A.vLen > 1;
        int mu_u = 0, mu_v = 0, mv_u = 0, mv_v = 0;
        bool mapping_ok = false;

        /*---------------------------------------------------------
          Phase 1: gather a few nearest-neighbor matches (hash lookup)
          until derive_uv_mapping succeeds or candidates are exhausted.
        ---------------------------------------------------------*/
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

            // Mapping inference requires multiple non-collinear samples.
            if (!mapping_ok && matches.size() >= 3) {
                mapping_ok = derive_uv_mapping(matches, need_u, need_v, mu_u, mu_v, mv_u, mv_v);
                if (mapping_ok)
                    break;
            }
        }

        /*---------------------------------------------------------
          Phase 2 (fast path): if (u,v) mapping is known, predict
          the corresponding (u,v) location in B for every A point
          and validate by distance.

          This avoids repeated hash searches, using FaceGrid::index
          as a direct mapping from (u,v) to FaceData index.

          FaceGrid::index stores FaceData indices (int) with size uLen*vLen.
        ---------------------------------------------------------*/
        if (mapping_ok) {
            const UVMatch& origin = matches.front();
            const FaceGrid& gridB = build.grids[static_cast<size_t>(faceB)];

            std::vector<MatchPair> local;
            local.reserve(static_cast<size_t>(countA));

            for (int ii = 0; ii < countA; ++ii) {
                int idx = build.face_cells[startA + static_cast<size_t>(ii)];

                // A relative offsets from mapping origin.
                int du = build.faces.u[idx] - origin.au;
                int dv = build.faces.v[idx] - origin.av;

                // Predicted B coordinates under linear mapping.
                int bu = origin.bu + du * mu_u + dv * mv_u;
                int bv = origin.bv + du * mu_v + dv * mv_v;

                // Reject out-of-range predicted coordinates.
                if (bu < 0 || bu >= B.uLen || bv < 0 || bv >= B.vLen)
                    continue;

                // Convert predicted (bu,bv) to the FaceData index for B.
                int b_idx = gridB.index[face_grid_index(gridB, bu, bv)];
                if (b_idx < 0)
                    continue;

                // Accept only if within tolerance.
                if (distance_sq(build.faces.x[idx], build.faces.y[idx], build.faces.z[idx],
                                build.faces.x[b_idx], build.faces.y[b_idx], build.faces.z[b_idx]) <= tol2)
                    local.push_back(MatchPair{idx, b_idx});
            }

            // If the mapping-based pass yielded any matches, use it and skip fallback.
            if (!local.empty()) {
                out.insert(out.end(), local.begin(), local.end());
                return;
            }
        }

        /*---------------------------------------------------------
          Phase 3 (fallback): hash lookup for every A point.

          This is used when mapping inference failed or mapping predicted
          no valid matches.
        ---------------------------------------------------------*/
        for (int ii = 0; ii < countA; ++ii) {
            int idx = build.face_cells[startA + static_cast<size_t>(ii)];
            int b_idx = find_match_in_face(hashB, build.faces,
                                           build.faces.x[idx], build.faces.y[idx], build.faces.z[idx],
                                           tol2, hash_cell);
            if (b_idx >= 0)
                out.push_back(MatchPair{idx, b_idx});
        }
    };

    /*-------------------------------------------------------------
      Process a contiguous span of candidate face pairs [begin,end),
      appending matches to 'out'.

      The implementation chooses the face with fewer points as the query side:
        - If faceA has more points than faceB, swap them for matching.
        - If swapped, swap indices back in the produced MatchPairs.
    -------------------------------------------------------------*/
    auto process_pairs = [&](size_t begin, size_t end, std::vector<MatchPair>& out) {
        for (size_t i = begin; i < end; ++i) {
            int faceA = pairs[i].a;
            int faceB = pairs[i].b;
            bool swapped = false;

            // Query the smaller face to reduce hash lookups / mapping work.
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

    /*-------------------------------------------------------------
      Thread selection:
        - threads<=0 -> use std::thread::hardware_concurrency()
        - enforce at least 1 thread
    -------------------------------------------------------------*/
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

        // Chunk candidate pairs across threads.
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

    /*-------------------------------------------------------------
      Merge per-thread match lists into a single vector and sort.
    -------------------------------------------------------------*/
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

/*=====================================================================
// resolve_matches_cgns
//
// Convert raw point matches into:
//   - boundary_records: unmatched face patches (candidate BC regions)
//   - conn_patches: matched face patches (candidate 1-to-1 connectivity)
//
// Inputs:
//   mesh     : provides Zone metadata and zone ordering (1-based IDs in output)
//   log      : prints “interior …” and “boundary …” summaries
//   build    : face grids + per-point metadata (block/face/i/j/k/u/v/x/y/z)
//   matches  : vector<MatchPair> linking FaceData indices that match
//
// Patch construction overview:
//
//   Step 1: Build NeighborInfo otherInfo[] indexed by FaceData index.
//           Each FaceData index gets at most one neighbor (first seen match).
//           otherInfo[idx].block == -1 => boundary point.
//
//   Step 2: For each zone boundary face:
//           Walk the FaceGrid in (u,v) to build maximal rectangular patches.
//
//           Patch growth rules:
//             - A patch grows in +u direction while neighbor info stays:
//                 * boundary across all those points, OR
//                 * same donor block+donor face and adjacent donor i/j/k
//             - Then grows in +v direction with the same constraints.
//
//           The “adjacent donor i/j/k” constraint requires that between
//           consecutive points in the growth direction exactly one of
//           |Δi|,|Δj|,|Δk| equals 1 (structured index step).
//
//   Step 3: For a completed patch:
//             - If boundary: emit BoundaryRecord(zone, face, begin/end, full?)
//             - If interior: compute CGNS transform and explicit donor ranges,
//               then emit ConnPatch.
//
// 2D handling:
//   - For zones with nk==1, recv and donor k are fixed to 1.
//   - transform[2] is forced to 0.
//   - The recv/donor point ranges are stored as IJK with k==1.
//=====================================================================*/
void resolve_matches_cgns(const Mesh& mesh, Logger& log,
                          FaceBuildResult& build,
                          const std::vector<MatchPair>& matches,
                          std::vector<BoundaryRecord>& boundary_records,
                          std::vector<ConnPatch>& conn_patches)
{
    // otherInfo.size() == number of FaceData points across all faces.
    // NeighborInfo default constructor (in connectivity_internal.hpp) sets block=-1.
    std::vector<NeighborInfo> otherInfo(build.faces.size());

    /*-------------------------------------------------------------
      Populate otherInfo from matches.

      The first match that references a FaceData index “claims” it.
      Subsequent matches touching that index are ignored.

      This restriction is required by the later patch grouping logic,
      which assumes each face point has either:
        - no neighbor (boundary)
        - exactly one neighbor (structured 1-to-1 interface)
    -------------------------------------------------------------*/
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

    /*-------------------------------------------------------------
      Iterate zones and faces, building rectangular patches from the
      per-point neighbor map.
    -------------------------------------------------------------*/
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

            // masterCount tracks how many (u,v) cells have been considered.
            // grid.consumed marks which cells have been assigned to a patch.
            int masterCount = 0;
            int uStart = 0;
            int vStart = 0;
            int uEnd = 0;
            int vEnd = 0;
            bool complete = false;

            while (!complete) {
                // Starting point of the current patch.                
                const auto startIdx = grid.index[face_grid_index(grid, uStart, vStart)];

                NeighborInfo curOther = otherInfo[startIdx];
                NeighborInfo startOther = curOther;

                // iDirIndex and jDirIndex capture the signed donor axis step
                // along receiver u and v directions respectively.
                // Values: ±1, ±2, ±3 representing ±I, ±J, ±K in donor.
                // Value 0 means “not yet inferred”.
                int iDirIndex = 0;
                int jDirIndex = 0;

                /*---------------------------------------------------------
                  Grow patch in +u direction.

                  Allowable continuation cases:
                    - both are boundary points (block==-1)
                    - both connect to same donor block+face and donor indices
                      advance by one step in one axis
                ---------------------------------------------------------*/
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

                    // Record the donor direction corresponding to receiver +u.
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

                /*---------------------------------------------------------
                  Grow patch in +v direction using the same constraints.
                ---------------------------------------------------------*/
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

                    // Record donor direction corresponding to receiver +v.
                    if (jDirIndex == 0 && nextOther.block != -1) {
                        if (delI == 1) jDirIndex = signed_axis(1, nextOther.i - curOther.i);
                        else if (delJ == 1) jDirIndex = signed_axis(2, nextOther.j - curOther.j);
                        else if (delK == 1) jDirIndex = signed_axis(3, nextOther.k - curOther.k);
                        else throw std::runtime_error("Unexpected J-direction delta");
                    }

                    curOther = nextOther;
                    vEnd = nextV;
                }

                /*---------------------------------------------------------
                  Mark all cells in the patch as consumed.
                ---------------------------------------------------------*/
                for (int uu = uStart; uu <= uEnd; ++uu)
                    for (int vv = vStart; vv <= vEnd; ++vv)
                        grid.consumed[face_grid_index(grid, uu, vv)] = 1;

                /*---------------------------------------------------------
                  Convert (u,v) patch bounds into CGNS/Plot3D point ranges.

                  fill_point_range() produces begin/end in 1-based IJK
                  and uses “+2” endpoints so that patch size in cells
                  maps to vertex ranges.

                  In 2D: k is pinned to 1 for both begin and end.
                ---------------------------------------------------------*/
                IJK recvBegin{};
                IJK recvEnd{};
                fill_point_range(z, face, uStart, vStart, uEnd, vEnd, recvBegin, recvEnd);

                if (z.nk() == 1) {
                    recvBegin[2] = 1;
                    recvEnd[2] = 1;
                }

                /*---------------------------------------------------------
                  Boundary patch: donor block == -1.
                ---------------------------------------------------------*/
                if (curOther.block == -1) {
                    bool isFull = (uStart == 0 && vStart == 0 && uEnd == uLen - 1 && vEnd == vLen - 1);
                    BoundaryRecord entry;
                    entry.zone = b + 1;                     // zones are 1-based in public output
                    entry.face = face_dir_from_id(face);
                    entry.full = isFull;
                    entry.begin = recvBegin;
                    entry.end = recvEnd;
                    boundary_records.push_back(entry);

                /*---------------------------------------------------------
                  Interior patch: has a donor block+face.

                  Steps:
                    - Determine missing iDirIndex/jDirIndex if only one was inferred
                      during growth (or neither was inferred for tiny patches).
                    - Build CGNS transform vector for receiver axes.
                    - Compute explicit donor point range using compute_donor_range().
                    - Emit ConnPatch.
                ---------------------------------------------------------*/
                } else {
                    int donorFace = curOther.face;

                    /*-----------------------------------------------------
                      If only one direction index was inferred, infer the other
                      based on donor face type (I/J/K face) so that two axes
                      span the donor face plane.

                      If neither was inferred (single-cell patches), choose
                      a default spanning set based on donorFace.
                    -----------------------------------------------------*/
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

                    /*-----------------------------------------------------
                      Build the 3-entry CGNS transform vector.

                      Convention:
                        - transform[axis_of_recv_i] = signed donor axis
                        - transform[axis_of_recv_j] = signed donor axis
                        - transform[axis_of_recv_k] = signed donor axis

                      Receiver face collapses one axis; set_normal_transform()
                      fills the collapsed-axis entry with the signed donor normal.
                    -----------------------------------------------------*/
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

                    /*-----------------------------------------------------
                      Compute donor begin/end point range.

                      donorStart is taken from the neighbor of the patch origin
                      in donor IJK space (1-based indices stored in FaceData).

                      compute_donor_range() applies the transform to the receiver
                      range extents and returns a consistent donor begin/end.

                      2D donor: pin k to 1 and force transform[2]=0 at write time.
                    -----------------------------------------------------*/
                    IJK donorBegin{};
                    IJK donorEnd{};
                    IJK donorStart = {startOther.i, startOther.j, startOther.k};
                    compute_donor_range(transform, recvBegin, recvEnd, donorStart, donorBegin, donorEnd);

                    const Zone& donorZ = mesh.zones()[curOther.block];
                    if (donorZ.nk() == 1) {
                        donorBegin[2] = 1;
                        donorEnd[2] = 1;
                    }

                    /*-----------------------------------------------------
                      Only write each interface once.

                      Ordering rule is (recvZone, recvFace) < (donorZone, donorFace).
                      This prevents duplicate reciprocal connectivity entries.
                    -----------------------------------------------------*/
                    int recvZone = b + 1;
                    int donorZone = curOther.block + 1;
                    if (should_write(recvZone, donorZone, face, donorFace)) {
                        int transform_cgns[3] = {transform[0], transform[1], transform[2]};

                        // If either side is 2D, CGNS transform's k component is 0.
                        if (z.nk() == 1 || donorZ.nk() == 1) {
                            transform_cgns[2] = 0;
                        }
                        bool isFull = (uStart == 0 && vStart == 0 && uEnd == uLen - 1 && vEnd == vLen - 1);
                        std::string recvFaceStr = facedir_to_string(face_dir_from_id(face));
                        std::string donorFaceStr = facedir_to_string(face_dir_from_id(donorFace));

                        // Log message includes explicit transform triplet.
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

                        // ConnPatch stores zone IDs as 1-based.
                        ConnPatch cp;
                        cp.recvZone = recvZone;
                        cp.donorZone = donorZone;
                        cp.recvDir = face_dir_from_id(face);
                        cp.donorDir = face_dir_from_id(donorFace);

                        // transform is stored as CGNS-style signed axes.
                        cp.transform[0] = transform_cgns[0];
                        cp.transform[1] = transform_cgns[1];
                        cp.transform[2] = transform_cgns[2];

                        // Store explicit receiver and donor ranges for writer.
                        cp.recvRange.begin = recvBegin;
                        cp.recvRange.end = recvEnd;
                        cp.donorRange.begin = donorBegin;
                        cp.donorRange.end = donorEnd;

                        conn_patches.push_back(cp);
                    }
                }

                /*---------------------------------------------------------
                  Choose next unconsumed (u,v) cell as the next patch origin.

                  masterCount increments through the flattened u-major order:
                    uu = (masterCount-1) % uLen
                    vv = (masterCount-1) / uLen

                  complete is set when all uLen*vLen cells are consumed.
                ---------------------------------------------------------*/
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

/*=====================================================================
// resolve_matches_plot3d
//
// Convert raw face-point matches on a Plot3D multi-block mesh into:
//   - boundary_records: unmatched face patches (candidate BC regions)
//   - conn_patches: matched face patches (candidate 1-to-1 connectivity)
//
// Inputs:
//   zones    : vector<Plot3DZone>, indexed by (block-1), providing ni/nj/nk
//              and implicit structured indexing for Plot3D blocks
//   log      : emits “interior …” and “boundary …” summaries
//   build    : FaceBuildResult produced by build_faces_plot3d():
//                * build.faces.*   : per-face-point metadata
//                * build.grids[]   : per-(block,face) FaceGrid (u×v layout)
//   matches  : vector<MatchPair> linking FaceData indices whose points
//              coincide geometrically within tolerance
//
// Outputs (appended to):
//   boundary_records : vector<BoundaryRecord>
//   conn_patches     : vector<ConnPatch>
//
// Patch construction overview:
//
//   Step 1: Build NeighborInfo otherInfo[] indexed by FaceData index.
//           Each FaceData index records at most one neighbor:
//             - block == -1 → no matching neighbor (boundary)
//             - otherwise   → donor block, donor face, and donor i/j/k
//           If a FaceData index appears in multiple MatchPairs, only the
//           first assignment is kept to enforce 1-to-1 connectivity.
//
//   Step 2: For each block and each of its 6 boundary faces:
//           Walk the FaceGrid (u,v) to form maximal rectangular patches.
//
//           Patch growth rules:
//             - Grow in +u direction while all points are compatible:
//                 * all boundary points, OR
//                 * all matched to the same donor block+face AND
//                   consecutive donor indices differ by exactly 1 in
//                   exactly one of (i,j,k)
//             - Then grow in +v direction with the same constraints.
//
//           During growth, the donor-space directions corresponding to
//           receiver U and V are inferred from donor i/j/k deltas.
//           These become iDirIndex and jDirIndex (signed axis indices).
//
//   Step 3: For a completed patch:
//             - Convert (uStart,vStart,uEnd,vEnd) into a receiver IJK
//               PointRange via fill_point_range().
//             - If all points are boundary:
//                   emit BoundaryRecord(zone, face, begin/end, full?)
//             - If all points match a donor:
//                   * finalize missing direction indices
//                   * build a CGNS-style transform[3]
//                   * compute donor PointRange via compute_donor_range()
//                   * emit ConnPatch with explicit recv/donor ranges
//
// Canonicalization:
//   should_write() enforces a consistent ordering so each interior
//   interface is emitted exactly once, even though matches are symmetric.
//
// 2D handling (nk == 1):
//   - recv and donor k indices are forced to 1
//   - transform[2] is forced to 0
//   - resulting PointRanges are valid 2D IJK ranges with k==1
//
// Notes:
//   - Plot3D output always uses explicit ranges and is written to text by bc_define_plot3d.cpp.
//=====================================================================*/
void resolve_matches_plot3d(const std::vector<Plot3DZone>& zones, Logger& log,
                            FaceBuildResult& build,
                            const std::vector<MatchPair>& matches,
                            std::vector<BoundaryRecord>& boundary_records,
                            std::vector<ConnPatch>& conn_patches)
{
    // Map each face-cell (FaceData entry) to its matched neighbor's metadata.
    //
    // otherInfo.size() == build.faces.size() (number of face-cells across all faces).
    // NeighborInfo fields:
    //   - block: 0-based block index (zone index in zones vector) or -1 for "no neighbor"
    //   - face:  0..5 face id on the donor block
    //   - i,j,k: 1-based vertex indices in donor block space for this matched face-cell
    //
    // Only a single neighbor is recorded per face-cell; subsequent matches involving
    // already-assigned cells are ignored to preserve a 1-to-1 pairing.
    std::vector<NeighborInfo> otherInfo(build.faces.size());

    // Populate otherInfo symmetrically from the flat list of MatchPair {a,b}.
    // match.a and match.b are indices into build.faces.* arrays.
    for (const auto& match : matches) {
        if (otherInfo[static_cast<size_t>(match.a)].block != -1 ||
            otherInfo[static_cast<size_t>(match.b)].block != -1)
            continue;

        // For face-cell "a", store the donor location derived from "b".
        otherInfo[match.a] = NeighborInfo{
            build.faces.block[match.b], build.faces.face[match.b],
            build.faces.i[match.b], build.faces.j[match.b], build.faces.k[match.b]};

        // For face-cell "b", store the donor location derived from "a".
        otherInfo[match.b] = NeighborInfo{
            build.faces.block[match.a], build.faces.face[match.a],
            build.faces.i[match.a], build.faces.j[match.a], build.faces.k[match.a]};
    }

    // Each Plot3D block produces up to 6 faces in FaceBuildResult, even when nk()==1.
    // The FaceBuildResult builder uses:
    //   face_id = block_index * 6 + face (face in [0..5])
    // The corresponding FaceGrid contains:
    //   - uLen, vLen: number of face-cells in each param direction (ints)
    //   - index:  size uLen*vLen, storing a FaceData index for each (u,v)
    //   - consumed: same size, marking which (u,v) have been processed into a patch
    const int nBlocks = static_cast<int>(zones.size());
    for (int b = 0; b < nBlocks; ++b)
    {
        const Plot3DZone& z = zones[b];
        for (int face = 0; face < 6; ++face)
        {
            int face_id = b * 6 + face;
            FaceGrid& grid = build.grids[static_cast<size_t>(face_id)];
            if (grid.index.empty()) continue;   // Face has no cells (degenerate or not present)

            int uLen = grid.uLen;   // face-cells in U (>=1)
            int vLen = grid.vLen;   // face-cells in V (==1 for 2D faces; >=1 for 3D)

            // masterCount is a linear scan counter over all uLen*vLen cells.
            // It advances until an unconsumed (u,v) is found, which becomes the next patch seed.
            int masterCount = 0;

            // Current patch seed and current patch extents in face-grid coordinates.
            // The extents [uStart..uEnd]x[vStart..vEnd] expand rightward (U) then upward (V)
            // while neighbor metadata remains compatible.
            int uStart = 0;
            int vStart = 0;
            int uEnd = 0;
            int vEnd = 0;

            bool complete = false;

            while (!complete) {
                // FaceData index at patch seed (uStart,vStart).
                const auto startIdx = grid.index[face_grid_index(grid, uStart, vStart)];

                // curOther tracks the neighbor at the current "end" of the growing patch.
                // startOther preserves the neighbor of the seed cell; it anchors donorStart
                // when computing the donor PointRange for the patch.
                NeighborInfo curOther = otherInfo[startIdx];
                NeighborInfo startOther = curOther;

                // iDirIndex/jDirIndex encode the donor-space axis that corresponds to the
                // receiver U and V directions across the face grid.
                //
                // Encoding:
                //   signed_axis(axis, delta) returns ±axis where axis is 1(I),2(J),3(K).
                // So iDirIndex and jDirIndex are in { -3,-2,-1,0,1,2,3 }.
                // 0 means "not yet determined".
                int iDirIndex = 0;
                int jDirIndex = 0;

                // ---------- Grow patch along U ----------
                // Attempt to extend the patch in U from uEnd to uEnd+1 (keeping vStart fixed).
                // Compatibility rules:
                //   - If both cells are boundary (block == -1), they are compatible.
                //   - If both cells have a neighbor, the neighbor block/face must match and
                //     the donor indices must be adjacent by exactly 1 in I or J or K.
                while (uEnd < uLen - 1) {
                    int nextU = uEnd + 1;
                    const auto idx = grid.index[face_grid_index(grid, nextU, vStart)];
                    NeighborInfo nextOther = otherInfo[idx];

                    int delI = std::abs(nextOther.i - curOther.i);
                    int delJ = std::abs(nextOther.j - curOther.j);
                    int delK = std::abs(nextOther.k - curOther.k);

                    bool ok = false;
                    if (nextOther.block == curOther.block && nextOther.block == -1) {
                        // boundary-adjacent: keep growing boundary patch
                        ok = true;
                    } else if (nextOther.block == curOther.block &&
                               nextOther.face == curOther.face &&
                               (delI == 1 || delJ == 1 || delK == 1)) {
                        // interior-adjacent within same donor face
                        ok = true;
                    }

                    if (!ok) break;

                    // Determine donor axis direction for receiver-U the first time an interior
                    // neighbor is encountered. The delta that equals 1 identifies which donor
                    // index changes along U, and the sign matches the direction.
                    if (iDirIndex == 0 && nextOther.block != -1) {
                        if (delI == 1) iDirIndex = signed_axis(1, nextOther.i - curOther.i);
                        else if (delJ == 1) iDirIndex = signed_axis(2, nextOther.j - curOther.j);
                        else if (delK == 1) iDirIndex = signed_axis(3, nextOther.k - curOther.k);
                        else throw std::runtime_error("Unexpected I-direction delta");
                    }

                    curOther = nextOther;
                    uEnd = nextU;
                }

                // Reset curOther to the neighbor at the current U edge, still on vStart,
                // before attempting to expand in V.
                curOther = otherInfo[grid.index[face_grid_index(grid, uEnd, vStart)]];

                // ---------- Grow patch along V ----------
                // Extend the patch in V from vEnd to vEnd+1 (keeping uEnd fixed).
                // The same compatibility rules apply as in U growth.
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

                    // Determine donor axis direction for receiver-V the first time an interior
                    // neighbor is encountered while expanding in V.
                    if (jDirIndex == 0 && nextOther.block != -1) {
                        if (delI == 1) jDirIndex = signed_axis(1, nextOther.i - curOther.i);
                        else if (delJ == 1) jDirIndex = signed_axis(2, nextOther.j - curOther.j);
                        else if (delK == 1) jDirIndex = signed_axis(3, nextOther.k - curOther.k);
                        else throw std::runtime_error("Unexpected J-direction delta");
                    }

                    curOther = nextOther;
                    vEnd = nextV;
                }

                // Mark every (u,v) cell inside the formed rectangle as consumed so it does
                // not seed another patch.
                //
                // consumed is stored as uint8_t in FaceGrid in this codebase
                // (1 byte per cell), and is written with literal 1/0.
                for (int uu = uStart; uu <= uEnd; ++uu)
                    for (int vv = vStart; vv <= vEnd; ++vv)
                        grid.consumed[face_grid_index(grid, uu, vv)] = 1;

                // Convert the face-grid rectangle into a Plot3D vertex PointRange.
                //
                // fill_point_range maps:
                //   - uStart/vStart/uEnd/vEnd (0-based cell indices)
                // to:
                //   - recvBegin/recvEnd (1-based vertex indices; inclusive bounds)
                //
                // The +2 in fill_point_range comes from converting cell extents to vertex
                // extents (cells count + 1 vertices).
                IJK recvBegin{};
                IJK recvEnd{};
                fill_point_range(z, face, uStart, vStart, uEnd, vEnd, recvBegin, recvEnd);

                // Plot3D 2D blocks have nk()==1. This code forces the K indices of the
                // resulting PointRange to 1 so downstream users treat it as 2D.
                if (z.nk() == 1) {
                    recvBegin[2] = 1;
                    recvEnd[2] = 1;
                }

                if (curOther.block == -1) {
                    // No neighbor found for this patch => boundary patch.
                    // full means the rectangle covers the entire face grid.
                    bool isFull = (uStart == 0 && vStart == 0 && uEnd == uLen - 1 && vEnd == vLen - 1);
                    BoundaryRecord entry;
                    entry.zone = b + 1;                     // convert 0-based block index to 1-based zone id
                    entry.face = face_dir_from_id(face);    // FaceDir enum derived from [0..5]
                    entry.full = isFull;
                    entry.begin = recvBegin;
                    entry.end = recvEnd;
                    boundary_records.push_back(entry);
                } else {
                    // Neighbor exists => interior 1-to-1 connectivity patch.

                    // donorFace is in [0..5] and is stored in NeighborInfo.
                    int donorFace = curOther.face;

                    // iDirIndex/jDirIndex can remain 0 if the patch is 1 cell wide in that
                    // direction (no adjacency step was observed). The following blocks
                    // fill missing direction(s) with a consistent axis choice based on the
                    // donor face type (I-face, J-face, or K-face).
                    //
                    // The choices are deterministic and use an axis not equal to the donor
                    // face's normal axis.
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
                        // Both are 0 (single-cell patch), or a defensive fallback.
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

                    // Build the CGNS-style 3-element transform vector for index mapping.
                    //
                    // Convention in this codebase:
                    //   transform[0] maps receiver I-axis
                    //   transform[1] maps receiver J-axis
                    //   transform[2] maps receiver K-axis
                    //
                    // For a receiver face, one of these axes is collapsed (constant) and
                    // corresponds to the face normal. The other two axes correspond to U and V.
                    int transform[3] = {0, 0, 0};
                    if (is_i_face(face)) {
                        transform[1] = iDirIndex;   // receiver J varies with U
                        transform[2] = jDirIndex;   // receiver K varies with V (or is collapsed in 2D)
                    } else if (is_j_face(face)) {
                        transform[0] = iDirIndex;   // receiver I varies with U
                        transform[2] = jDirIndex;   // receiver K varies with V
                    } else {
                        transform[0] = iDirIndex;   // receiver I varies with U
                        transform[1] = jDirIndex;   // receiver J varies with V
                    }

                    // Fill the remaining transform component associated with the receiver
                    // face normal. The sign depends on whether receiver and donor are on
                    // opposite sides (+/-) and on CGNS 1-to-1 transform conventions.
                    set_normal_transform(transform, face, donorFace);

                    // Compute donor range endpoints from:
                    //   - receiver range (recvBegin..recvEnd)
                    //   - donor seed index (startOther.i/j/k)
                    //   - transform vector
                    //
                    // compute_donor_range returns donorBegin/donorEnd as IJK indices in
                    // the donor zone (still 1-based).
                    IJK donorBegin{};
                    IJK donorEnd{};
                    IJK donorStart = {startOther.i, startOther.j, startOther.k};
                    compute_donor_range(transform, recvBegin, recvEnd, donorStart, donorBegin, donorEnd);

                    // 2D Plot3D blocks use K==1; enforce that invariant for donor as well.
                    const Plot3DZone& donorZ = zones[static_cast<size_t>(curOther.block)];
                    if (donorZ.nk() == 1) {
                        donorBegin[2] = 1;
                        donorEnd[2] = 1;
                    }

                    // Convert indices into the exported ConnPatch format:
                    //   - recvZone/donorZone are 1-based
                    //   - recvDir/donorDir are FaceDir enums
                    int recvZone = b + 1;
                    int donorZone = curOther.block + 1;

                    // should_write enforces a canonical ordering so each interface is
                    // emitted once (not twice from both sides).
                    if (should_write(recvZone, donorZone, face, donorFace)) {
                        // For Plot3D outputs, the same integer transform convention is used
                        // as the CGNS writer in this codebase.
                        int transform_cgns[3] = {transform[0], transform[1], transform[2]};

                        // For 2D blocks, CGNS-style transform uses 0 for the collapsed axis.
                        if (z.nk() == 1 || donorZ.nk() == 1) {
                            transform_cgns[2] = 0;
                        }
                        bool isFull = (uStart == 0 && vStart == 0 && uEnd == uLen - 1 && vEnd == vLen - 1);
                        std::string recvFaceStr = facedir_to_string(face_dir_from_id(face));
                        std::string donorFaceStr = facedir_to_string(face_dir_from_id(donorFace));

                        // Logging is split to distinguish full-face connections from
                        // partial patches; both include the integer transform.
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

                        // ConnPatch stores explicit vertex ranges (recvRange/donorRange).
                        // For Plot3D output, the explicit ranges are the authoritative representation
                        // and are written by bc_define_plot3d.cpp to the "1to1s" text file.
                        ConnPatch cp;
                        cp.recvZone = recvZone;
                        cp.donorZone = donorZone;
                        cp.recvDir = face_dir_from_id(face);
                        cp.donorDir = face_dir_from_id(donorFace);

                        // Store transform and explicit ranges.
                        cp.transform[0] = transform_cgns[0];
                        cp.transform[1] = transform_cgns[1];
                        cp.transform[2] = transform_cgns[2];

                        cp.recvRange.begin = recvBegin;
                        cp.recvRange.end = recvEnd;
                        cp.donorRange.begin = donorBegin;
                        cp.donorRange.end = donorEnd;

                        conn_patches.push_back(cp);
                    }
                }

                // Advance masterCount until a not-yet-consumed cell is found, then use it
                // as the next patch seed. Linear indexing:
                //   uu = (masterCount-1) % uLen
                //   vv = (masterCount-1) / uLen
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
                
                // When every cell is consumed, finish this face and move to the next.
                if (masterCount == uLen * vLen) {
                    complete = true;
                }
            }
        }
    }
}

} // namespace fs::conn
