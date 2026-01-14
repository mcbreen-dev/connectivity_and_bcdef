
#ifndef CONNECTIVITY_HPP
#define CONNECTIVITY_HPP

#include "common.hpp"
#include "connectivity_writer.hpp"
#include <vector>

namespace fs {
    class Mesh;
    struct Plot3DZone;
    class Logger;

    struct BoundaryPatch {
        int     zone = 0;           // 1-based zone index
        FaceDir face = FaceDir::IMIN;
        IJK     vtxBegin{};
        IJK     vtxEnd{};
        bool    fullFace = false;
    };

    struct ConnectivityResult {
        std::vector<BoundaryPatch> boundary;
        std::vector<ConnPatch> connections;
    };

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
