//─────────────────────────────────────────────────────────────
// File: src/connectivity_hash.cpp
// Spatial hashing and candidate generation
//─────────────────────────────────────────────────────────────
#include "connectivity_internal.hpp"

#include <algorithm>

namespace fs::conn {

HashBuildResult build_hash_candidates(const FaceBuildResult& build,
                                      const std::vector<FaceSummary>& summaries,
                                      double tol)
{
    HashBuildResult result;
    size_t face_count = summaries.size();
    result.hash_cell = (tol > 0.0) ? tol : 1.0;
    result.face_hashes.resize(face_count);

    for (size_t face_id = 0; face_id < face_count; ++face_id) {
        auto& hash = result.face_hashes[face_id].cells;
        size_t start = build.cell_start[face_id];
        int count = build.cell_count[face_id];
        hash.reserve(static_cast<size_t>(count) * 2 + 1);
        for (int ii = 0; ii < count; ++ii) {
            int idx = build.face_cells[start + static_cast<size_t>(ii)];
            hash[cell_key(build.faces.x[idx], build.faces.y[idx], build.faces.z[idx], result.hash_cell)]
                .push_back(idx);
        }
    }

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
    double bucket = std::max(tol * 4.0, extent / 32.0);
    if (bucket <= 0.0)
        bucket = 1.0;

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
