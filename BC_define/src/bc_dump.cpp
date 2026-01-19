//─────────────────────────────────────────────────────────────
// File: src/bc_dump.cpp
// CGNS zone / connectivity / BC dump utility
//─────────────────────────────────────────────────────────────
//
// This file implements a small command-line program that *reads* a CGNS
// mesh and prints a human-readable summary of:
//
//   - Zone dimensions (Ni, Nj, Nk)
//   - Existing GridConnectivity1to1_t entries (“1to1” abutting interfaces)
//   - Existing BC_t entries (boundary conditions)
//
// The intent is diagnostic / inspection:
//   - Verify that bc_define (and/or other tools) wrote the expected
//     connectivity and BCs.
//   - Quickly inspect 2-D vs 3-D handling (Nk==1 indicates 2-D zones).
//
// Output conventions:
//   - Zones are printed as "Zone1", "Zone2", ... for readability,
//     independent of the raw CGNS Zone_t names.
//   - When parsing connectivity names, we expect the writer’s naming
//     pattern:  Z<recv>_<recvDir>_to_Z<donor>_<donorDir>
//     If tokens are missing, we fall back to the donor zone name stored
//     by CGNS, and/or to Z<idx>.
//   - Ranges are printed in a consistent 3-axis (I,J,K) format.
//     For 2-D zones (Nk==1) we pin K=1.
//
// Notes / assumptions:
//   - This tool uses the CGNS C API directly; error checking is best-effort
//     since it is a dump tool (it skips malformed entries where possible).
//   - It does not modify the CGNS file.
//─────────────────────────────────────────────────────────────
#include "mesh_io.hpp"
#include "mesh_utils.hpp"
#include "logger.hpp"
#include "common.hpp"

#include <iomanip>
#include <iostream>
#include <string>
#include <unordered_map>
#include <vector>

extern "C" {
    #include <cgnslib.h>
}


/*=====================================================================
  Global helpers

  These are used by the dump routines to:
    - translate raw CGNS zone names into stable short labels ("Z<idx>")
    - determine 2-D vs 3-D behavior per zone via Nk (Nk==1 => 2-D)

  They are initialized in main() after opening the mesh.
=====================================================================*/
static std::unordered_map<std::string,std::string> g_zname2nice;
static std::vector<long long> g_zoneNk;   // Nk per zone (index-1)

#ifndef CGNS_MAX_NAME_LENGTH
#define CGNS_MAX_NAME_LENGTH 32   // fallback for Homebrew CGNS
#endif


/*=====================================================================
  bc_type_string

  Provide a concise label for a BCType_t for printing.

  The dump focuses on a few common BC types. Everything else is grouped
  into "other" to keep output compact.
=====================================================================*/
static std::string bc_type_string(CGNS_ENUMT(BCType_t) t)
{
    switch (t) {
        case CGNS_ENUMV(BCWallViscousHeatFlux): return "wall-viscous";
        case CGNS_ENUMV(BCWallInviscid):        return "wall-inviscid";
        case CGNS_ENUMV(BCInflow):              return "inflow";
        case CGNS_ENUMV(BCOutflow):             return "outflow";
        case CGNS_ENUMV(BCSymmetryPlane):       return "symmetry";
        case CGNS_ENUMV(BCFarfield):            return "farfield";
        default:                                return "other";
    }
}

/*=====================================================================
  nice_zone_name

  Convert a raw CGNS Zone_t name into a stable "Zone<N>" label.

  Currently this is unconditional: it always returns "Zone<idx>".
  (Keeping the raw name can be useful in some workflows; if desired,
  this could be changed to only rename zones that match a pattern.)
=====================================================================*/
static std::string nice_zone_name(const std::string& raw, int idx)
{
    // Always print Zone1, Zone2, … regardless of the raw CGNS name
    // (you can make this conditional on a pattern if you prefer)
    return "Zone" + std::to_string(idx);
}

/*=====================================================================
  dump_1to1

  Print all GridConnectivity1to1_t entries for a single Zone.

  Inputs:
    fn, B : CGNS file + base IDs
    Z     : 1-based zone index in CGNS (matches Zone.idx)

  Behavior:
    - Uses cg_n1to1 / cg_1to1_read to read all 1-to-1 links.
    - Attempts to decode the link name into receiver/donor zone + face
      strings based on the naming convention used by connectivity_writer.
    - Prints receiver point range (normalized into 3D printing format)
      and the integer transform vector T{t1,t2,t3}.

  2D handling:
    - If Nk==1 for the receiver zone, the range returned by CGNS is
      treated as [imin, jmin, imax, jmax] (2D), and K is set to 1 for
      printing.
=====================================================================*/
static void dump_1to1(int fn, int B, int Z)
{
    int nLink = 0;
    cg_n1to1(fn, B, Z, &nLink);
    if (nLink == 0) return;

    std::cout << "    Interior 1to1 connections (" << nLink << "):\n";

    for (int i = 1; i <= nLink; ++i)
    {
        /* CGNS buffers for reading one GridConnectivity1to1_t entry. */
        char      name  [CGNS_MAX_NAME_LENGTH + 1] = {};
        char      donor [CGNS_MAX_NAME_LENGTH + 1] = {};
        cgsize_t  range[6]        = {0};   /* receiver  begin / end */
        cgsize_t  donor_range[6]  = {0};   /* donor     begin / end */
        int       T[3]            = {0};   /* transform */

        if (cg_1to1_read(fn, B, Z, i,
                         name, donor,
                         range, donor_range, T) != CG_OK)
            continue;

        /*----------------------------------------------------------
          Decode writer-style name tokens.

          Expected name:
            Z<recv>_<recvDir>_to_Z<donor>_<donorDir>

          The dump tool is robust to deviations: if decoding fails, it
          falls back to:
            - receiver zone index Z
            - donor zone name from CGNS (or g_zname2nice mapping)
        ----------------------------------------------------------*/
        std::string sname(name);
        std::vector<std::string> tok;
        size_t start = 0, pos;
        while ((pos = sname.find('_', start)) != std::string::npos) {
            tok.push_back(sname.substr(start, pos - start));
            start = pos + 1;
        }
        tok.push_back(sname.substr(start));

        std::string recvZ   = tok.size() > 0 ? tok[0] : ("Z" + std::to_string(Z));
        std::string recvDir = tok.size() > 1 ? tok[1] : "";
        std::string donorZ  = (tok.size() > 3 ? tok[3]
                                              : (g_zname2nice.count(donor)
                                                    ? g_zname2nice[donor]
                                                    : donor));
        std::string donorDir = tok.size() > 4 ? tok[4] : "";

        /*----------------------------------------------------------
          Normalize receiver range for printing.

          In CGNS 2D paths, the 1-to-1 API returns a 4-value range
          (imin,jmin,imax,jmax) in range[0..3]. We print as 3D by
          pinning K=1.
        ----------------------------------------------------------*/
        bool twoD = (g_zoneNk[Z-1] == 1);      // Nk==1 → 2-D zone

        cgsize_t rmin[3], rmax[3];
        if (twoD) {
            rmin[0] = range[0];   rmin[1] = range[1];   rmin[2] = 1;
            rmax[0] = range[2];   rmax[1] = range[3];   rmax[2] = 1;
        } else {
            rmin[0] = range[0];   rmin[1] = range[1];   rmin[2] = range[2];
            rmax[0] = range[3];   rmax[1] = range[4];   rmax[2] = range[5];
        }

        /* Print a single summary line per link. */
        std::cout << "      • Donor " << donorZ << ' ' << donorDir
                    << " → Receiver "   << recvZ  << ' ' << recvDir
                    << "   r:[" << rmin[0] << ',' << rmin[1] << ',' << rmin[2]
                    << "]-["  << rmax[0] << ',' << rmax[1] << ',' << rmax[2] << ']'
                    << "  T{" << T[0]   << ',' << T[1]   << ',' << T[2]   << "}\n";
    }
}


/*=====================================================================
  dump_bcs

  Print all BC_t entries for a single Zone.

  Inputs:
    fn, B : CGNS file + base IDs
    Z     : 1-based zone index

  Behavior:
    - Uses cg_nbocos / cg_boco_info / cg_boco_read to enumerate BCs.
    - If a BC uses FamilySpecified, attempts to read its FamilyName via
      cg_goto + cg_famname_read.
    - If a BC is a PointRange with exactly two points, prints its range.

  2D handling:
    - If Nk==1, ranges are interpreted as [imin,jmin,imax,jmax] and K is
      pinned to 1 for printing.
=====================================================================*/
static void dump_bcs(int fn,int B,int Z)
{
    int nBC = 0;
    cg_nbocos(fn,B,Z,&nBC);
    if (nBC == 0) return;

    std::cout << "    BCs (" << nBC << "):\n";

    /* Cache BC metadata so we can do a clean second pass for printing. */
    struct BCEntry {
        std::string name;
        std::string family;
        CGNS_ENUMT(BCType_t) type = CGNS_ENUMV(BCTypeNull);
        CGNS_ENUMT(PointSetType_t) pset = CGNS_ENUMV(PointSetTypeNull);
        cgsize_t nPnts = 0;
        cgsize_t range[6] = {0};
        bool has_range = false;
    };

    std::vector<BCEntry> entries;
    entries.reserve(static_cast<size_t>(nBC));

    /* First pass: collect entries. */
    for (int ibc=1; ibc<=nBC; ++ibc)
    {
        char name[CGNS_MAX_NAME_LENGTH+1] = {};
        CGNS_ENUMT(BCType_t)        type;
        CGNS_ENUMT(PointSetType_t)  pset;
        cgsize_t nPnts = 0;
        int NormalInd[3] = {0,0,0};
        cgsize_t NormalListSz = 0;
        CGNS_ENUMT(DataType_t) NormalDataType;
        int nDataset = 0;

        cg_boco_info(fn,B,Z,ibc,name,&type,&pset,&nPnts,
                      NormalInd,&NormalListSz,&NormalDataType,&nDataset);

        BCEntry entry;
        entry.name = name;
        entry.type = type;
        entry.pset = pset;
        entry.nPnts = nPnts;

        /* Read the PointRange (if applicable). */
        if (pset == CGNS_ENUMV(PointRange) && nPnts == 2) {
            cg_boco_read(fn,B,Z,ibc,entry.range,nullptr);
            entry.has_range = true;
        }

        /* If FamilySpecified, try to read the FamilyName string. */
        char famname[CGNS_MAX_NAME_LENGTH+1] = {};
        if (cg_goto(fn, B, "Zone_t", Z, "ZoneBC_t", 1, "BC_t", ibc, "end") == CG_OK) {
            if (cg_famname_read(famname) == CG_OK)
                entry.family = famname;
        }

        entries.push_back(entry);
    }

    /* Second pass: pretty-print. */
    for (int ibc=1; ibc<=nBC; ++ibc)
    {
        const BCEntry& entry = entries[static_cast<size_t>(ibc - 1)];
        std::string label = entry.family.empty() ? entry.name : entry.family;

        std::cout << "      • " << label;

        /* If BC is not FamilySpecified, print its explicit CGNS type. */
        if (entry.type != CGNS_ENUMV(FamilySpecified))
            std::cout << "  (" << bc_type_string(entry.type) << ")";

        /* If range exists, normalize to 3D-style output and print it. */
        if (entry.has_range)
        {
            bool twoD = (g_zoneNk[Z-1] == 1);
            cgsize_t rmin[3] = {0, 0, 0};
            cgsize_t rmax[3] = {0, 0, 0};
            if (twoD) {
                rmin[0] = entry.range[0]; rmin[1] = entry.range[1]; rmin[2] = 1;
                rmax[0] = entry.range[2]; rmax[1] = entry.range[3]; rmax[2] = 1;
            } else {
                rmin[0] = entry.range[0]; rmin[1] = entry.range[1]; rmin[2] = entry.range[2];
                rmax[0] = entry.range[3]; rmax[1] = entry.range[4]; rmax[2] = entry.range[5];
            }
            std::cout << "  ["<<rmin[0]<<','<<rmin[1]<<','<<rmin[2]<<"]-["
                      << rmax[0]<<','<<rmax[1]<<','<<rmax[2]<<']';
        }
        std::cout << '\n';
    }
}


/*=====================================================================
  main

  Usage:
      bc_dump mesh.cgns

  Flow:
    1) Open CGNS mesh read-only
    2) Build global helpers:
         - g_zname2nice: raw zone name -> "Z<idx>"
         - g_zoneNk:     Nk per zone (for 2D detection)
    3) Print zone summary table
    4) For each zone:
         - dump_1to1() : print 1-to-1 connectivity links
         - dump_bcs()  : print BC entries
    5) Close mesh and exit

  Exit codes:
    0 : success
    1 : CLI usage error
    2 : runtime error (exception)
=====================================================================*/
int main(int argc,char** argv)
{
    if (argc!=2){ std::cerr << "Usage: bc_dump mesh.cgns\n"; return 1; }

    const std::string path = argv[1];
    fs::Logger log("bc_dump.log");
    log.info("Opening "+path);

    try{
        fs::Mesh mesh; mesh.open(path,false);

        /* Map native CGNS zone names to "Z<idx>" labels for cleaner printing. */
        g_zname2nice.clear();
        for (const auto& z : mesh.zones()) {
            g_zname2nice[z.name] = "Z" + std::to_string(z.idx);
        }

        /* Record Nk per zone so dump routines can detect 2-D zones. */
        g_zoneNk.clear();
        for (const auto& z : mesh.zones())
            g_zoneNk.push_back(z.nk());
            
        /* Print a compact zone dimension table. */
        std::cout << "Zones in file \""<<path<<"\":\n";
        std::cout << "Name      Ni   Nj   Nk\n"
                     "-------- ---- ---- ----\n";

        for (const auto &z : mesh.zones()){
            std::string nice = "Zone" + std::to_string(z.idx);          // always ≤8 chars
            //std::cout   << std::setw(4) << z.idx       << ' '          // Idx
            std::cout   << std::left  << std::setw(8)  << nice         // Name
                        << std::left << ' '
                        << std::setw(4) << z.ni() << ' '               // Ni
                        << std::setw(4) << z.nj() << ' '               // Nj
                        << std::setw(4) << z.nk() << '\n';             // Nk

            dump_1to1(mesh.file_id(),mesh.base_id(),z.idx);
            dump_bcs (mesh.file_id(),mesh.base_id(),z.idx);
        }

        mesh.close();
        log.info("Finished dump.");
        return 0;
    }
    catch(const std::exception& ex){
        log.error(ex.what());
        std::cerr << ex.what() << '\n';
        return 2;
    }
}
