/*
  File: include/mesh_utils.hpp

  CGNS tree cleanup utilities used by bc_define when --overwrite is enabled.

  This header declares two functions implemented in src/mesh_utils.cpp:

    - purge_BC_and_connectivity:
        Deletes, for every Zone_t, the entire "ZoneBC" and "ZoneGridConnectivity"
        child nodes (if present). This removes all BC_t and 1to1 connectivity
        information in one operation per node.

        Called in:
          - bc_define_main.cpp when overwrite==true before running connectivity/BC write.

    - prune_unused_families:
        Scans Zone_t and BC_t nodes for referenced FamilyName_t strings and removes
        any Family_t nodes in the base that are not referenced.

        Called in:
          - bc_define_main.cpp when overwrite==true after writing new BCs.

  Type / size notes:
    - Mesh holds the CGNS file handle (int) and base id (int); these are used to
      navigate the CGNS tree with cg_goto / cg_delete_node.
    - Logger is used for emitting status lines; Logger itself holds an ofstream
      and a mutex (see logger.hpp).
*/
#pragma once

#include "mesh_io.hpp"
#include "logger.hpp"

extern "C" {
    #include <cgnslib.h>
}

namespace fs {

    /*
      purge_BC_and_connectivity(mesh, log)

      Effects (src/mesh_utils.cpp):
        - For each zone index Z in mesh.zones():
            cg_goto(fn, B, "Zone_t", Z, NULL)
            cg_delete_node("ZoneBC")
            cg_delete_node("ZoneGridConnectivity")
        - cg_delete_node return codes treated as:
            0               : deleted successfully
            CG_NODE_NOT_FOUND: node absent (no-op)
            CG_ERROR        : treated as absent in this code (covers common "H5Gopen failed")
            otherwise       : prints CGNS error and throws std::runtime_error

      Logging:
        - Writes one INFO line: "Purged existing BC / connectivity nodes"

      Called in:
        - bc_define_main.cpp when --overwrite flag is passed
    */
    void purge_BC_and_connectivity(Mesh& mesh, Logger& log);

    /*
      prune_unused_families(mesh, log)

      Effects (src/mesh_utils.cpp):
        - Builds a set of referenced family names by:
            1) checking Zone_t for a FamilyName_t (cg_famname_read at "Zone_t", Z)
            2) checking each BC_t for a FamilyName_t (cg_famname_read under "ZoneBC_t"/"BC_t")
        - Enumerates Family_t nodes in the base (cg_nfamilies / cg_family_read)
        - Deletes any Family_t nodes whose names are not in the referenced set:
            cg_goto(fn, B, "end")
            cg_delete_node(family_name)
        - Logs each deleted family as:
            "Removed unused Family_t <name>"

      Called in:
        - bc_define_main.cpp when --overwrite flag is passed (after writing BCs)
    */
    void prune_unused_families(Mesh& mesh, Logger& log);
}
