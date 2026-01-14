//─────────────────────────────────────────────────────────────
// File: src/connectivity_core.cpp
// ConnectivityDetector entry points
//─────────────────────────────────────────────────────────────
#include "connectivity.hpp"
#include "connectivity_internal.hpp"
#include "connectivity_writer.hpp"
#include "logger.hpp"

namespace fs {

ConnectivityResult ConnectivityDetector::run_collect(
    Mesh& mesh, Logger& log, double tol, bool show_progress, int threads)
{
    std::vector<conn::BoundaryRecord> boundary_records;
    std::vector<ConnPatch> conn_patches;

    const int fn = mesh.file_id();
    const int B  = mesh.base_id();

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

    size_t totalFaces = 0;
    if (show_progress) {
        for (const auto& z : mesh.zones())
            totalFaces += (z.nk() > 1 ? 6 : 4);
    }

    conn::FaceBuildResult build = conn::build_faces_cgns(mesh, show_progress);

    if (show_progress) {
        log.info("Faces read: " + std::to_string(build.faces_loaded) + "/" +
                 std::to_string(totalFaces));
    }

    if (show_progress) {
        log.info("Creating 1:1 connections");
    }

    if (build.faces.empty()) {
        log.info("Connectivity scan complete.");
        return {};
    }

    std::vector<conn::FaceSummary> summaries = conn::build_face_summaries(build);
    conn::HashBuildResult hash = conn::build_hash_candidates(build, summaries, tol);

    boundary_records.reserve(64);
    std::vector<conn::MatchPair> matches =
        conn::collect_match_pairs(build, summaries, hash, tol, threads);
    conn::resolve_matches_cgns(mesh, log, build, matches,
                               boundary_records, conn_patches);

    log.info("Connectivity scan complete.");
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

ConnectivityResult ConnectivityDetector::run_plot3d(
    const std::vector<Plot3DZone>& zones, Logger& log, double tol, int threads)
{
    std::vector<conn::BoundaryRecord> boundary_records;
    std::vector<ConnPatch> conn_patches;

    const int nBlocks = static_cast<int>(zones.size());
    if (nBlocks == 0) {
        return {};
    }

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

    size_t totalFaces = 0;
    for (const auto& z : zones)
        totalFaces += (z.nk() > 1 ? 6 : 4);

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

std::vector<BoundaryPatch> ConnectivityDetector::run(
    Mesh& mesh, Logger& log, double tol, bool show_progress, int threads)
{
    ConnectivityResult result = run_collect(mesh, log, tol, show_progress, threads);
    for (const auto& cp : result.connections) {
        write_1to1(mesh, cp);
    }
    return result.boundary;
}

} // namespace fs
