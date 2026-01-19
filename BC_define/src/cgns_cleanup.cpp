/*─────────────────────────────────────────────────────────────
  File: src/cgns_cleanup.cpp

  CGNS mesh maintenance utilities used by bc_define --overwrite.

  This file provides two related cleanup operations:

    1) purge_BC_and_connectivity(mesh, log)
       - Removes existing boundary condition and connectivity containers:
           * ZoneBC_t               (named "ZoneBC" in CGNS files)
           * ZoneGridConnectivity_t (named "ZoneGridConnectivity")
       - Implemented by deleting the *entire* child node under each Zone_t,
         which is simpler and faster than iterating and deleting each BC_t
         or GridConnectivity_t entry individually.

    2) prune_unused_families(mesh, log)
       - After overwrite, removes Family_t nodes under the Base that are no
         longer referenced by any Zone_t or BC_t node.
       - This keeps the CGNS file tidy and avoids stale family definitions.

  Notes:
    - These routines assume "base_id() == 1" (first base), consistent with
      Mesh::open in mesh_io.cpp.
    - The CGNS library API is used directly. Errors are checked via CG_CALL
      (see common.hpp) where appropriate, and explicit CGNS error codes are
      handled where deletion semantics vary by backend (HDF5).
─────────────────────────────────────────────────────────────*/
#include "cgns_cleanup.hpp"
#include "common.hpp"
#include <unordered_set>
#include <vector>

extern "C" {
    #include <cgnslib.h>
}

#ifndef CGNS_MAX_NAME_LENGTH
#define CGNS_MAX_NAME_LENGTH 32
#endif

namespace bcdef {

/*=====================================================================
  delete_child

  Delete a named child node under a given Zone_t node.

  Parameters:
    fn, B : CGNS file and base ids
    Z     : 1-based zone index
    childName : name of the child node to delete (e.g., "ZoneBC")

  Behavior:
    - Positions the CGNS cursor at Zone_t(Z)
    - Attempts to delete the named child via cg_delete_node(childName)

  Deletion semantics / tolerated error codes:
    - ier == 0              : deleted successfully
    - CG_NODE_NOT_FOUND     : node does not exist; treated as success
    - CG_ERROR              : some HDF5 backends report a failed open
                             (H5Gopen) as CG_ERROR; interpreted here as
                             "child absent" and treated as success

  Any other return code is considered a real error and throws.
=====================================================================*/
static void delete_child(int fn, int B, int Z, const char* childName)
{
    /* position at Zone_t */
    CG_CALL(cg_goto(fn, B, "Zone_t", Z, NULL),
            "CGNS: cg_goto Zone_t failed");

    /* attempt to delete the child node */
    int ier = cg_delete_node(childName);

    if (ier == 0)                       return;     // deleted OK
    if (ier == CG_NODE_NOT_FOUND)       return;     // child absent → fine
    if (ier == CG_ERROR)             return;        // H5Gopen failed → child absent
    
    /* anything else means a real problem */
    cg_error_print();
    throw std::runtime_error("CGNS: cg_delete_node failed (" +
                                 std::string(childName) + ")");

}

/*=====================================================================
  purge_BC_and_connectivity

  Overwrite helper: remove existing BC/connectivity containers from all zones.

  This deletes, per Zone_t:
    - ZoneBC (ZoneBC_t)
    - ZoneGridConnectivity (ZoneGridConnectivity_t)

  Deleting these containers removes all their children (BC_t, GridConnectivity_t,
  GridConnectivity1to1_t, etc.) in one operation.
=====================================================================*/
void purge_BC_and_connectivity(Mesh& mesh, Logger& log)
{
    int fn = mesh.file_id(), B = mesh.base_id();

    for (const auto& z : mesh.zones()) {
        int Z = z.idx;
        // Pointwise & CGNS defaults
        delete_child(fn, B, Z, "ZoneBC");
        delete_child(fn, B, Z, "ZoneGridConnectivity");
    }
    
    log.info("Purged existing BC / connectivity nodes");
}

/*=====================================================================
  add_family_if_present

  If BC_t(bc_index) under Zone_t(Z) has a FamilyName, add it to 'families'.

  This is used to determine which Family_t nodes are still referenced
  after an overwrite cleanup.

  Notes:
    - cg_goto positions at .../ZoneBC_t/BC_t(bc_index)
    - cg_famname_read reads the FamilyName_t value for the current node
=====================================================================*/
static void add_family_if_present(std::unordered_set<std::string>& families,
                                  int fn, int B, int Z, int bc_index)
{
    char famname[CGNS_MAX_NAME_LENGTH + 1] = {};
    if (cg_goto(fn, B, "Zone_t", Z, "ZoneBC_t", 1, "BC_t", bc_index, "end") != CG_OK)
        return;
    if (cg_famname_read(famname) == CG_OK)
        families.insert(famname);
}

/*=====================================================================
  add_zone_family_if_present

  If Zone_t(Z) itself has a FamilyName, add it to 'families'.

  Some CGNS files tag zones with a FamilyName_t, independent of BCs.
=====================================================================*/
static void add_zone_family_if_present(std::unordered_set<std::string>& families,
                                       int fn, int B, int Z)
{
    char famname[CGNS_MAX_NAME_LENGTH + 1] = {};
    if (cg_goto(fn, B, "Zone_t", Z, "end") != CG_OK)
        return;
    if (cg_famname_read(famname) == CG_OK)
        families.insert(famname);
}

/*=====================================================================
  prune_unused_families

  After overwrite, remove Family_t nodes that are no longer referenced.

  A Family_t is considered "referenced" if its name appears as:
    - Zone_t/FamilyName_t
    - ZoneBC_t/BC_t/FamilyName_t

  Implementation steps:
    1) Build a set of referenced family names by scanning all zones and BCs
    2) Read all Family_t nodes under the Base
    3) Delete families not present in the referenced set

  Notes / limitations:
    - This only checks direct FamilyName links from zones and BCs.
      If your CGNS file references families from other node types,
      they are not considered here.
=====================================================================*/
void prune_unused_families(Mesh& mesh, Logger& log)
{
    int fn = mesh.file_id(), B = mesh.base_id();

    /*---------------------------------------------------------
      Step 1: Collect referenced family names
    ---------------------------------------------------------*/
    std::unordered_set<std::string> referenced;

    for (const auto& z : mesh.zones()) {
        add_zone_family_if_present(referenced, fn, B, z.idx);

        int nBC = 0;
        if (cg_nbocos(fn, B, z.idx, &nBC) != CG_OK)
            continue;
        for (int bc = 1; bc <= nBC; ++bc)
            add_family_if_present(referenced, fn, B, z.idx, bc);
    }

    /*---------------------------------------------------------
      Step 2: Enumerate families declared under the base
    ---------------------------------------------------------*/
    int nFamilies = 0;
    if (cg_nfamilies(fn, B, &nFamilies) != CG_OK)
        return;

    std::vector<std::string> to_delete;
    to_delete.reserve(static_cast<size_t>(nFamilies));

    char name[CGNS_MAX_NAME_LENGTH + 1];
    int nboco = 0, ngeos = 0;
    for (int f = 1; f <= nFamilies; ++f) {
        if (cg_family_read(fn, B, f, name, &nboco, &ngeos) != CG_OK)
            continue;
        if (referenced.find(name) == referenced.end())
            to_delete.emplace_back(name);
    }
    
    /*---------------------------------------------------------
      Step 3: Delete unreferenced families
    ---------------------------------------------------------*/
    for (const auto& fname : to_delete) {
        if (cg_goto(fn, B, "end") != CG_OK)
            continue;
        if (cg_delete_node(fname.c_str()) == CG_OK)
            log.info("Removed unused Family_t " + fname);
    }
}

}
