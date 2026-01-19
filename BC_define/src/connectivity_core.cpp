//─────────────────────────────────────────────────────────────
// File: src/connectivity_core.cpp
//
// ConnectivityDetector entry points.
//
// This file provides the public-facing implementations of
// ConnectivityDetector::run*, which orchestrate the full
// connectivity detection pipeline for both CGNS and Plot3D meshes.
//
// Responsibilities in this file:
//   - Validate mesh handedness
//   - Build face representations for all zones
//   - Drive face matching and connectivity resolution
//   - Collect boundary patches and 1-to-1 connectivity patches
//   - Optionally write CGNS GridConnectivity1to1 nodes
//
// All geometric operations, hashing, and matching logic live in
// connectivity_face_build.cpp, connectivity_hash.cpp,
// connectivity_verify.cpp and connectivity_utils.cpp.
//─────────────────────────────────────────────────────────────
#include "connectivity.hpp"
#include "connectivity_internal.hpp"
#include "connectivity_writer.hpp"
#include "logger.hpp"

namespace bcdef {

/*=====================================================================
  ConnectivityDetector::run_collect  (CGNS)

  Full connectivity detection for a CGNS mesh, without writing any
  connectivity nodes back to the file.

  This function:
    - Detects and validates zone handedness
    - Builds face geometry for all zones
    - Matches coincident faces across zones
    - Resolves matches into:
        * ConnPatch entries (1-to-1 connections)
        * BoundaryRecord entries (unmatched boundary faces)

  The returned ConnectivityResult contains:
    - boundary   : faces that require BC assignment
    - connections: face pairs suitable for GridConnectivity1to1_t

  Writing of CGNS connectivity nodes is handled by run().
=====================================================================*/
ConnectivityResult ConnectivityDetector::run_collect(
    Mesh& mesh, Logger& log, double tol, bool show_progress, int threads)
{
    /*---------------------------------------------------------
      Temporary storage used during matching and resolution
    ---------------------------------------------------------*/
    std::vector<conn::BoundaryRecord> boundary_records;
    std::vector<ConnPatch> conn_patches;

    const int fn = mesh.file_id();
    const int B  = mesh.base_id();

    /*---------------------------------------------------------
      Handedness validation:
      All 3D zones must be right-handed for consistent face
      orientation and transform computation.

      2D zones (nk == 1) are excluded from this check.
    ---------------------------------------------------------*/
    bool checked_handedness = false;
    for (const auto& z : mesh.zones()) {
        if (z.nk() == 1) continue;
        checked_handedness = true;
        if (conn::is_left_handed(fn, B, z))
            throw std::runtime_error("Zone " + std::to_string(z.idx) +
                                     " left-handed");
    }

    if (checked_handedness) {
        log.info("All zones right-handed (OK)");
    }

    /*---------------------------------------------------------
      Progress reporting: total number of logical faces
      (6 for 3D zones, 4 for 2D zones).
    ---------------------------------------------------------*/
    size_t totalFaces = 0;
    if (show_progress) {
        for (const auto& z : mesh.zones())
            totalFaces += (z.nk() > 1 ? 6 : 4);
    }

    /*---------------------------------------------------------
      Build face geometry for all zones.

      This reads coordinate data from CGNS and constructs
      canonical face representations suitable for hashing
      and comparison.
    ---------------------------------------------------------*/
    conn::FaceBuildResult build = conn::build_faces_cgns(mesh, show_progress);

    if (show_progress) {
        log.info("Faces read: " + std::to_string(build.faces_loaded) + "/" +
                 std::to_string(totalFaces));
    }

    if (show_progress) {
        log.info("Creating 1:1 connections");
    }

    /*---------------------------------------------------------
      Early exit if no faces were constructed.
    ---------------------------------------------------------*/
    if (build.faces.empty()) {
        log.info("Connectivity scan complete.");
        return {};
    }

    /*---------------------------------------------------------
      Reduce face data into compact summaries used for hashing
      and candidate matching.
    ---------------------------------------------------------*/
    std::vector<conn::FaceSummary> summaries = conn::build_face_summaries(build);

    /*---------------------------------------------------------
      Build spatial hash buckets of candidate face pairs based
      on bounding boxes and tolerance.
    ---------------------------------------------------------*/
    conn::HashBuildResult hash = conn::build_hash_candidates(build, summaries, tol);

    /*---------------------------------------------------------
      Perform detailed face matching in parallel (if enabled),
      producing matched face pairs.
    ---------------------------------------------------------*/
    boundary_records.reserve(64);
    std::vector<conn::MatchPair> matches =
        conn::collect_match_pairs(build, summaries, hash, tol, threads);
    
    /*---------------------------------------------------------
      Resolve matches:
        - Generate ConnPatch entries for matched faces
        - Record unmatched faces as BoundaryRecord entries
    ---------------------------------------------------------*/
    conn::resolve_matches_cgns(mesh, log, build, matches,
                               boundary_records, conn_patches);

    log.info("Connectivity scan complete.");

    /*---------------------------------------------------------
      Log boundary faces that will later receive BCs.
    ---------------------------------------------------------*/
    if (!boundary_records.empty()) {
        log.info("Assigning boundary conditions");
        for (const auto& entry : boundary_records) {
            std::string recvFaceStr = facedir_to_string(entry.face);
            if (entry.full) {
                log.info("boundary Z" + std::to_string(entry.zone) + " " + recvFaceStr);
            } else {
                log.info(
                    "boundary PATCH Z" + std::to_string(entry.zone) + " " + recvFaceStr +
                    "  r:[" + std::to_string(entry.begin[0]) + "," +
                    std::to_string(entry.begin[1]) + "," +
                    std::to_string(entry.begin[2]) + "]-[" +
                    std::to_string(entry.end[0]) + "," +
                    std::to_string(entry.end[1]) + "," +
                    std::to_string(entry.end[2]) + "]"
                );
            }
        }
    }

    /*---------------------------------------------------------
      Convert internal BoundaryRecord entries into public
      BoundaryPatch structures.
    ---------------------------------------------------------*/
    std::vector<BoundaryPatch> boundary_patches;
    boundary_patches.reserve(boundary_records.size());
    for (const auto& entry : boundary_records) {
        BoundaryPatch bp;
        bp.zone = entry.zone;
        bp.face = entry.face;
        bp.vtxBegin = entry.begin;
        bp.vtxEnd = entry.end;
        bp.fullFace = entry.full;
        boundary_patches.push_back(bp);
    }

    ConnectivityResult result;
    result.boundary = std::move(boundary_patches);
    result.connections = std::move(conn_patches);
    return result;
}

/*=====================================================================
  ConnectivityDetector::run_plot3d

  Connectivity detection for Plot3D multi-block meshes.

  Differences from CGNS path:
    - Uses Plot3DZone instead of Mesh/Zone
    - No file I/O; results are written as text files elsewhere
    - No overwrite or progress toggles
=====================================================================*/
ConnectivityResult ConnectivityDetector::run_plot3d(
    const std::vector<Plot3DZone>& zones, Logger& log, double tol, int threads)
{
    std::vector<conn::BoundaryRecord> boundary_records;
    std::vector<ConnPatch> conn_patches;

    const int nBlocks = static_cast<int>(zones.size());
    if (nBlocks == 0) {
        return {};
    }

    /*---------------------------------------------------------
      Handedness validation for Plot3D blocks.
    ---------------------------------------------------------*/
    bool checked_handedness = false;
    for (const auto& z : zones) {
        if (z.nk() == 1) continue;
        checked_handedness = true;
        if (conn::is_left_handed_plot3d(z))
            throw std::runtime_error("Zone " + std::to_string(z.idx) +
                                     " left-handed");
    }

    if (checked_handedness) {
        log.info("All zones right-handed (OK)");
    }

    /*---------------------------------------------------------
      Count logical faces for logging.
    ---------------------------------------------------------*/
    size_t totalFaces = 0;
    for (const auto& z : zones)
        totalFaces += (z.nk() > 1 ? 6 : 4);

    /*---------------------------------------------------------
      Build face geometry from Plot3D coordinate arrays.
    ---------------------------------------------------------*/
    conn::FaceBuildResult build = conn::build_faces_plot3d(zones);

    log.info("Faces read: " + std::to_string(build.faces_loaded) + "/" +
             std::to_string(totalFaces));
    log.info("Creating 1:1 connections");

    if (build.faces.empty()) {
        log.info("Connectivity scan complete.");
        return {};
    }

    std::vector<conn::FaceSummary> summaries = conn::build_face_summaries(build);
    conn::HashBuildResult hash = conn::build_hash_candidates(build, summaries, tol);

    boundary_records.reserve(64);
    std::vector<conn::MatchPair> matches =
        conn::collect_match_pairs(build, summaries, hash, tol, threads);
    conn::resolve_matches_plot3d(zones, log, build, matches,
                                 boundary_records, conn_patches);

    log.info("Connectivity scan complete.");

    /*---------------------------------------------------------
      Log boundary patches.
    ---------------------------------------------------------*/
    if (!boundary_records.empty()) {
        log.info("Assigning boundary conditions");
        for (const auto& entry : boundary_records) {
            std::string recvFaceStr = facedir_to_string(entry.face);
            if (entry.full) {
                log.info("boundary Z" + std::to_string(entry.zone) + " " + recvFaceStr);
            } else {
                log.info(
                    "boundary PATCH Z" + std::to_string(entry.zone) + " " + recvFaceStr +
                    "  r:[" + std::to_string(entry.begin[0]) + "," +
                    std::to_string(entry.begin[1]) + "," +
                    std::to_string(entry.begin[2]) + "]-[" +
                    std::to_string(entry.end[0]) + "," +
                    std::to_string(entry.end[1]) + "," +
                    std::to_string(entry.end[2]) + "]"
                );
            }
        }
    }

    std::vector<BoundaryPatch> boundary_patches;
    boundary_patches.reserve(boundary_records.size());
    for (const auto& entry : boundary_records) {
        BoundaryPatch bp;
        bp.zone = entry.zone;
        bp.face = entry.face;
        bp.vtxBegin = entry.begin;
        bp.vtxEnd = entry.end;
        bp.fullFace = entry.full;
        boundary_patches.push_back(bp);
    }

    ConnectivityResult result;
    result.boundary = std::move(boundary_patches);
    result.connections = std::move(conn_patches);
    return result;
}

/*=====================================================================
  ConnectivityDetector::run

  Convenience wrapper for CGNS meshes.

  This function:
    - Runs full connectivity detection via run_collect
    - Writes GridConnectivity1to1_t nodes to the CGNS file
    - Returns only boundary patches for BC assignment
=====================================================================*/
std::vector<BoundaryPatch> ConnectivityDetector::run(
    Mesh& mesh, Logger& log, double tol, bool show_progress, int threads)
{
    ConnectivityResult result = run_collect(mesh, log, tol, show_progress, threads);
    
    /*---------------------------------------------------------
      Write CGNS 1-to-1 connectivity entries.
    ---------------------------------------------------------*/
    for (const auto& cp : result.connections) {
        write_1to1(mesh, cp);
    }
    return result.boundary;
}

} // namespace bcdef
