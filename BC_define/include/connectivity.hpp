#ifndef CONNECTIVITY_HPP
#define CONNECTIVITY_HPP

/*
  File: include/connectivity.hpp

  Public interface for the structured connectivity detector.

  Usage in this codebase:

    - bc_define_main.cpp (CGNS path):
        * fs::ConnectivityDetector::run(mesh, log, fs::GEOM_TOL, overwrite, threads)
          returns boundary patches and also writes 1to1 connections into the CGNS
          file (via write_1to1 for each connection).

    - bc_define_main.cpp (Plot3D path):
        * fs::ConnectivityDetector::run_plot3d(zones, log, fs::GEOM_TOL, threads)
          returns:
            - boundary patches (used to generate Plot3D BC ranges)
            - connection patches (written to a text "1to1s" file)

    - bc_define_cgns.cpp:
        * write_boundary_conditions(...) consumes BoundaryPatch data:
            zone / face / vtxBegin / vtxEnd / fullFace

  Size / type notes (typical LP64 platforms):
    - int: 32-bit signed
    - long long: 64-bit signed
    - FaceDir: enum class backed by uint8_t (1 byte) in common.hpp
    - IJK: std::array<long long, 3> => 3 * 8 = 24 bytes
    - BoundaryPatch: contains:
        int zone (4) + FaceDir (1) + padding + two IJK (2*24) + bool (1) + padding
      Actual sizeof depends on alignment, but payload is dominated by the two IJK arrays.
*/

#include "common.hpp"
#include "connectivity_writer.hpp"
#include <vector>

namespace fs {
    class Mesh;
    struct Plot3DZone;
    class Logger;

    /*
      BoundaryPatch

      Represents a boundary region on a zone face that is NOT matched to another
      zone face as a 1:1 interface.

      Produced by:
        - ConnectivityDetector::run_collect(...) (CGNS)
        - ConnectivityDetector::run_plot3d(...)  (Plot3D)

      Consumed by:
        - bc_define_cgns.cpp:
            write_boundary_conditions(mesh, patches, ...)
          which converts vtxBegin/vtxEnd (vertex indexing) into a FaceCenter
          PointRange for BC writing.

        - bc_define_plot3d.cpp:
            build_plot3d_bcs(...) converts vtxBegin/vtxEnd into FaceCenter
            ranges for writing a Plot3D "bcs" file.

      Field semantics:
        - zone:
            1-based zone index (matches fs::Zone::idx in mesh_io.hpp).
            Type: int (CGNS zone numbering is int-based).

        - face:
            Which face of the zone: IMIN/IMAX/JMIN/JMAX/KMIN/KMAX.
            Type: FaceDir (uint8_t underlying).

        - vtxBegin / vtxEnd:
            Vertex-index endpoints defining the patch region on the face.
            Type: IJK = std::array<long long,3>.
            Indices are 1-based in the same convention as CGNS point indices.
            For 2D zones (Nk == 1), K is kept as 1 by the connectivity code.

            Note: These are VERTEX-space indices. BC writing uses FACE-CENTER
            indices, so bc_define_helpers.cpp provides face_center_range(...)
            which maps vertex indices to FaceCenter ranges by subtracting 1 from
            the varying axes.

        - fullFace:
            True when the patch covers the entire valid face-cell grid:
              (uStart==0 && vStart==0 && uEnd==uLen-1 && vEnd==vLen-1)
            as detected in connectivity_verify.cpp.
            This is used for logging and for downstream logic that may treat
            full-face patches differently.
    */
    struct BoundaryPatch {
        int     zone = 0;           // 1-based zone index
        FaceDir face = FaceDir::IMIN;
        IJK     vtxBegin{};
        IJK     vtxEnd{};
        bool    fullFace = false;
    };

    /*
      ConnectivityResult

      Bundle return type used by:
        - ConnectivityDetector::run_collect(...) (CGNS)
        - ConnectivityDetector::run_plot3d(...)  (Plot3D)

      Fields:
        - boundary:
            List of BoundaryPatch records (unmatched regions).

        - connections:
            List of ConnPatch records describing matched 1:1 interfaces.

          CGNS path:
            - Produced by resolve_matches_cgns(...)
            - Written into the CGNS file by ConnectivityDetector::run(...)
              using write_1to1(mesh, cp) for each entry.

          Plot3D path:
            - Produced by resolve_matches_plot3d(...)
            - Not written into a CGNS file; instead bc_define_plot3d.cpp writes
              a text "1to1s" file using the explicit ranges stored in ConnPatch.
    */
    struct ConnectivityResult {
        std::vector<BoundaryPatch> boundary;
        std::vector<ConnPatch> connections;
    };

    /*
      ConnectivityDetector

      Top-level entry points that coordinate the connectivity scan.

      Implementation files involved:
        - connectivity_core.cpp:
            Implements run, run_collect, run_plot3d.
        - connectivity_face_build.cpp:
            Reads face geometry and builds FaceBuildResult (cell-center points).
        - connectivity_hash.cpp:
            Builds spatial hash and candidate face pairs.
        - connectivity_verify.cpp:
            Confirms matches and constructs BoundaryRecord and ConnPatch outputs.
        - connectivity_writer.cpp:
            Writes CGNS GridConnectivity1to1_t when run(...) is used.

      Methods:

        run(mesh, logger, tol, show_progress, threads) -> vector<BoundaryPatch>
          - Calls run_collect(...) to detect boundary + connections.
          - Writes each connection into the CGNS file via write_1to1(mesh, cp).
          - Returns only the boundary patches (used by bc_define to write BCs).
          - show_progress affects:
              * handedness scan logging
              * face count logging
              * "Creating 1:1 connections" logging
            (see connectivity_core.cpp)

        run_collect(mesh, logger, tol, show_progress, threads) -> ConnectivityResult
          - CGNS-only detection path.
          - Validates handedness for 3D zones by sampling a local Jacobian
            determinant (connectivity_utils.cpp::is_left_handed).
          - Builds face point clouds (connectivity_face_build.cpp::build_faces_cgns).
          - Builds summaries and candidate pairs, then matches and resolves:
              build_face_summaries(...)
              build_hash_candidates(...)
              collect_match_pairs(...)
              resolve_matches_cgns(...)
          - Converts BoundaryRecord -> BoundaryPatch and returns both outputs.

        run_plot3d(zones, logger, tol, threads) -> ConnectivityResult
          - Plot3D-only detection path.
          - Same algorithmic steps as run_collect, but geometry is read from
            Plot3DZone coordinate arrays.
          - Returned connections are used for text outputs, not CGNS writes.
    */
    class ConnectivityDetector {
    public:
        static std::vector<BoundaryPatch> run(
            Mesh& mesh, Logger& logger, double tol, bool show_progress = true,
            int threads = 0);
        static ConnectivityResult run_collect(
            Mesh& mesh, Logger& logger, double tol, bool show_progress = true,
            int threads = 0);
        static ConnectivityResult run_plot3d(
            const std::vector<Plot3DZone>& zones, Logger& logger, double tol,
            int threads = 0);
    };
}

#endif // CONNECTIVITY_HPP
