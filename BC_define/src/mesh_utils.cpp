/*─────────────────────────────────────────────────────────────
  mesh_utils.cpp  –  bulk-delete BC & connectivity for --overwrite
  Strategy: remove entire ZoneBC_t and ZoneGridConnectivity_t nodes
            instead of iterating through their children.
─────────────────────────────────────────────────────────────*/
#include "mesh_utils.hpp"
#include "common.hpp"
#include <unordered_set>
#include <vector>

extern "C" {
    #include <cgnslib.h>
}

#ifndef CGNS_MAX_NAME_LENGTH
#define CGNS_MAX_NAME_LENGTH 32
#endif

namespace fs {

/* small helper: delete child if it exists */
static void delete_child(int fn, int B, int Z, const char* childName)
{
    /* position at Zone_t */
    CG_CALL(cg_goto(fn, B, "Zone_t", Z, NULL),
            "CGNS: cg_goto Zone_t failed");

    /* attempt to delete the child node */
    int ier = cg_delete_node(childName);

    if (ier == 0)                       return;     // deleted OK
    if (ier == CG_NODE_NOT_FOUND)       return;     // child absent → fine
    if (ier == CG_ERROR)             return;     // H5Gopen failed → child absent
    
    /* anything else means a real problem */
    cg_error_print();
    throw std::runtime_error("CGNS: cg_delete_node failed (" +
                                 std::string(childName) + ")");

}

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

static void add_family_if_present(std::unordered_set<std::string>& families,
                                  int fn, int B, int Z, int bc_index)
{
    char famname[CGNS_MAX_NAME_LENGTH + 1] = {};
    if (cg_goto(fn, B, "Zone_t", Z, "ZoneBC_t", 1, "BC_t", bc_index, "end") != CG_OK)
        return;
    if (cg_famname_read(famname) == CG_OK)
        families.insert(famname);
}

static void add_zone_family_if_present(std::unordered_set<std::string>& families,
                                       int fn, int B, int Z)
{
    char famname[CGNS_MAX_NAME_LENGTH + 1] = {};
    if (cg_goto(fn, B, "Zone_t", Z, "end") != CG_OK)
        return;
    if (cg_famname_read(famname) == CG_OK)
        families.insert(famname);
}

void prune_unused_families(Mesh& mesh, Logger& log)
{
    int fn = mesh.file_id(), B = mesh.base_id();

    std::unordered_set<std::string> referenced;

    for (const auto& z : mesh.zones()) {
        add_zone_family_if_present(referenced, fn, B, z.idx);

        int nBC = 0;
        if (cg_nbocos(fn, B, z.idx, &nBC) != CG_OK)
            continue;
        for (int bc = 1; bc <= nBC; ++bc)
            add_family_if_present(referenced, fn, B, z.idx, bc);
    }

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

    for (const auto& fname : to_delete) {
        if (cg_goto(fn, B, "end") != CG_OK)
            continue;
        if (cg_delete_node(fname.c_str()) == CG_OK)
            log.info("Removed unused Family_t " + fname);
    }
}

}
