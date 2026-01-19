//─────────────────────────────────────────────────────────────
// File: src/connectivity_hash.cpp
// Spatial hashing and candidate face-pair generation
//─────────────────────────────────────────────────────────────
//
// This file builds a coarse spatial acceleration structure used to
// reduce the number of face–face comparisons required during
// connectivity detection.
//
// The hashing is deliberately conservative: it may generate false
// positives (face pairs that do not actually match), but it must not
// miss any true matches. All candidates produced here are later
// verified geometrically and topologically in collect_match_pairs().

#include "connectivity_internal.hpp"

#include <algorithm>

namespace fs::conn {

/*=====================================================================
// build_hash_candidates
//
// Build spatial hash structures and generate candidate face pairs that
// may be geometrically adjacent and therefore eligible for 1-to-1
// connectivity.
//
// Inputs:
//   build     : FaceBuildResult containing per-face grids and per-point
//               FaceData indices (x,y,z stored in build.faces)
//   summaries : one FaceSummary per face, containing bounding boxes,
//               centroid, dimensions, and face identifiers
//   tol       : geometric tolerance used both for hashing resolution
//               and later bounds overlap checks
//
// Outputs (in HashBuildResult):
//   hash_cell   : spatial cell size used for point hashing
//   face_hashes : per-face hash tables mapping CellKey → FaceData indices
//                 (used later during point-level matching)
//   pairs       : vector<FacePair> (a,b) giving candidate face IDs whose
//                 bounding boxes overlap within tolerance
//
// Algorithm overview:
//
//   Phase 1: Per-face point hashing
//     - For each face, iterate over all FaceData indices belonging to
//       that face.
//     - Each face point is mapped to a CellKey computed by quantizing
//       its (x,y,z) coordinates using hash_cell.
//     - The resulting hash (CellKey → list of FaceData indices) is
//       stored per face. This is used later to efficiently find nearby
//       points across faces.
//
//   Phase 2: Face bounding-box bucketing
//     - Compute a global bounding box over all faces.
//     - Choose a bucket size based on mesh extent and tolerance.
//     - Each face bounding box is quantized into a 3D range of bucket
//       indices and inserted into all buckets it overlaps.
//     - Each bucket stores a list of face IDs whose bounds intersect it.
//
//   Phase 3: Candidate face-pair generation
//     - For each bucket, generate all unique (i,j) face pairs inside
//       that bucket.
//     - Deduplicate face pairs globally.
//     - Perform a final bounds-overlap test with tolerance.
//     - Store surviving pairs in result.pairs.
//
// Notes on tolerances and robustness:
//   - hash_cell is set to tol if tol > 0, otherwise 1.0.
//   - Bucket size is max(4*tol, global_extent/32), ensuring that faces
//     that are close but lie on opposite sides of a bucket boundary are
//     still co-bucketed.
//   - No assumptions are made about face orientation or normal; this
//     stage is purely spatial.
//
// Performance characteristics:
//   - Hashing is O(Npoints) for face points.
//   - Bucketing is O(Nfaces × buckets-per-face), typically small.
//   - Candidate generation avoids the O(Nfaces²) explosion for large
//     multi-block meshes.
//
//=====================================================================*/
HashBuildResult build_hash_candidates(const FaceBuildResult& build,
                                      const std::vector<FaceSummary>& summaries,
                                      double tol)
{
    HashBuildResult result;

    size_t face_count = summaries.size();

    // Spatial quantization size for point hashing.
    // If tol <= 0, fall back to unit-sized cells.
    result.hash_cell = (tol > 0.0) ? tol : 1.0;

    // One hash table per face: CellKey → vector<FaceData index>
    result.face_hashes.resize(face_count);

    // ------------------------------------------------------------------
    // Phase 1: Build per-face point hash tables
    // ------------------------------------------------------------------
    for (size_t face_id = 0; face_id < face_count; ++face_id) {
        auto& hash = result.face_hashes[face_id].cells;

        size_t start = build.cell_start[face_id];
        int count = build.cell_count[face_id];

        // Reserve roughly 2× entries to reduce rehashing
        hash.reserve(static_cast<size_t>(count) * 2 + 1);
        for (int ii = 0; ii < count; ++ii) {
            int idx = build.face_cells[start + static_cast<size_t>(ii)];
            hash[cell_key(build.faces.x[idx], build.faces.y[idx], build.faces.z[idx], result.hash_cell)]
                .push_back(idx);
        }
    }

    // ------------------------------------------------------------------
    // Phase 2: Compute global bounds and face buckets
    // ------------------------------------------------------------------
    FaceBounds global_bounds;
    for (const auto& summary : summaries) {
        if (!summary.bounds.valid)
            continue;
        expand_bounds(global_bounds, summary.bounds.min);
        expand_bounds(global_bounds, summary.bounds.max);
    }

    double extent_x = global_bounds.max.x - global_bounds.min.x;
    double extent_y = global_bounds.max.y - global_bounds.min.y;
    double extent_z = global_bounds.max.z - global_bounds.min.z;
    double extent = std::max(extent_x, std::max(extent_y, extent_z));

    // Bucket size used for face-level spatial hashing
    double bucket = std::max(tol * 4.0, extent / 32.0);
    if (bucket <= 0.0)
        bucket = 1.0;

    // CellKey → list of face IDs whose bounding boxes overlap the cell
    std::unordered_map<CellKey, std::vector<int>, CellKeyHash> face_buckets;
    face_buckets.reserve(face_count * 2);

    for (const auto& summary : summaries) {
        if (!summary.bounds.valid)
            continue;
        const FaceBounds& bnd = summary.bounds;
        std::int64_t ix0 = quantize_coord_offset(bnd.min.x, global_bounds.min.x, bucket);
        std::int64_t iy0 = quantize_coord_offset(bnd.min.y, global_bounds.min.y, bucket);
        std::int64_t iz0 = quantize_coord_offset(bnd.min.z, global_bounds.min.z, bucket);
        std::int64_t ix1 = quantize_coord_offset(bnd.max.x, global_bounds.min.x, bucket);
        std::int64_t iy1 = quantize_coord_offset(bnd.max.y, global_bounds.min.y, bucket);
        std::int64_t iz1 = quantize_coord_offset(bnd.max.z, global_bounds.min.z, bucket);

        for (std::int64_t ix = ix0; ix <= ix1; ++ix)
            for (std::int64_t iy = iy0; iy <= iy1; ++iy)
                for (std::int64_t iz = iz0; iz <= iz1; ++iz)
                    face_buckets[CellKey{ix, iy, iz}].push_back(summary.face_id);
    }

    // ------------------------------------------------------------------
    // Phase 3: Generate and filter candidate face pairs
    // ------------------------------------------------------------------
    std::vector<FacePair> candidates;
    for (const auto& entry : face_buckets) {
        const auto& bucket_faces = entry.second;
        for (size_t i = 0; i < bucket_faces.size(); ++i) {
            for (size_t j = i + 1; j < bucket_faces.size(); ++j) {
                int a = bucket_faces[i];
                int b = bucket_faces[j];
                if (a == b)
                    continue;
                if (a > b)
                    std::swap(a, b);
                candidates.push_back(FacePair{a, b});
            }
        }
    }

    // Sort and deduplicate candidate pairs
    std::sort(candidates.begin(), candidates.end(),
              [](const FacePair& lhs, const FacePair& rhs) {
                  if (lhs.a != rhs.a) return lhs.a < rhs.a;
                  return lhs.b < rhs.b;
              });
    candidates.erase(std::unique(candidates.begin(), candidates.end(),
                                 [](const FacePair& lhs, const FacePair& rhs) {
                                     return lhs.a == rhs.a && lhs.b == rhs.b;
                                 }),
                     candidates.end());
                     
    // Final bounds-overlap filter with tolerance
    result.pairs.reserve(candidates.size());
    for (const auto& pair : candidates) {
        const FaceSummary& a = summaries[static_cast<size_t>(pair.a)];
        const FaceSummary& b = summaries[static_cast<size_t>(pair.b)];
        if (!bounds_overlap(a.bounds, b.bounds, tol))
            continue;
        result.pairs.push_back(pair);
    }

    return result;
}

} // namespace fs::conn
