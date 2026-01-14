#pragma once
#include "mesh_io.hpp"
#include "logger.hpp"

extern "C" {
    #include <cgnslib.h>
}

namespace fs {
    void purge_BC_and_connectivity(Mesh& mesh, Logger& log);
    void prune_unused_families(Mesh& mesh, Logger& log);
}
