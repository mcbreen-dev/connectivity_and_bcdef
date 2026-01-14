/*─────────────────────────────────────────────────────────────
  File: src/bc_dump.cpp  — CGNS 4.5 zone / connectivity / BC dump
─────────────────────────────────────────────────────────────*/
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

/* global lookup so helper functions can see readable names */
static std::unordered_map<std::string,std::string> g_zname2nice;
static std::vector<long long> g_zoneNk;   // Nk per zone (index-1)

#ifndef CGNS_MAX_NAME_LENGTH
#define CGNS_MAX_NAME_LENGTH 32   // fallback for Homebrew CGNS
#endif

/*------------------------------------------------------------*/
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

/*------------------------------------------------------------*/
static std::string nice_zone_name(const std::string& raw, int idx)
{
    // Always print Zone1, Zone2, … regardless of the raw CGNS name
    // (you can make this conditional on a pattern if you prefer)
    return "Zone" + std::to_string(idx);
}

/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
static void dump_1to1(int fn, int B, int Z)
{
    int nLink = 0;
    cg_n1to1(fn, B, Z, &nLink);
    if (nLink == 0) return;

    std::cout << "    Interior 1to1 connections (" << nLink << "):\n";

    for (int i = 1; i <= nLink; ++i)
    {
        /* ── CGNS call buffers ───────────────────────────── */
        char      name  [CGNS_MAX_NAME_LENGTH + 1] = {};
        char      donor [CGNS_MAX_NAME_LENGTH + 1] = {};
        cgsize_t  range[6]        = {0};   /* receiver  begin / end */
        cgsize_t  donor_range[6]  = {0};   /* donor     begin / end */
        int       T[3]            = {0};   /* transform */

        if (cg_1to1_read(fn, B, Z, i,
                         name, donor,
                         range, donor_range, T) != CG_OK)
            continue;

        /* -------------------------------------------------- *
         *  Decode name:  Z<recv>_<recvDir>_to_Z<don>_<donDir>
         * -------------------------------------------------- */
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

        /* -------------------------------------------------- *
         *  Print: always show the 6-value receiver range     *
         * -------------------------------------------------- */

        /* ---------- build nice receiver range for printing ---------- */
        bool twoD = (g_zoneNk[Z-1] == 1);      // Nk==1 → 2-D zone

        cgsize_t rmin[3], rmax[3];
        if (twoD) {
            /* mapping: range = [imin jmin imax jmax 0 0] */
            rmin[0] = range[0];   rmin[1] = range[1];   rmin[2] = 1;
            rmax[0] = range[2];   rmax[1] = range[3];   rmax[2] = 1;
        } else {
            rmin[0] = range[0];   rmin[1] = range[1];   rmin[2] = range[2];
            rmax[0] = range[3];   rmax[1] = range[4];   rmax[2] = range[5];
        }

        std::cout << "      • Donor " << donorZ << ' ' << donorDir
                    << " → Receiver "   << recvZ  << ' ' << recvDir
                    << "   r:[" << rmin[0] << ',' << rmin[1] << ',' << rmin[2]
                    << "]-["  << rmax[0] << ',' << rmax[1] << ',' << rmax[2] << ']'
                    << "  T{" << T[0]   << ',' << T[1]   << ',' << T[2]   << "}\n";
    }
}

/*------------------------------------------------------------*/
static void dump_bcs(int fn,int B,int Z)
{
    int nBC = 0;
    cg_nbocos(fn,B,Z,&nBC);
    if (nBC == 0) return;

    std::cout << "    BCs (" << nBC << "):\n";

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

        if (pset == CGNS_ENUMV(PointRange) && nPnts == 2) {
            cg_boco_read(fn,B,Z,ibc,entry.range,nullptr);
            entry.has_range = true;
        }

        char famname[CGNS_MAX_NAME_LENGTH+1] = {};
        if (cg_goto(fn, B, "Zone_t", Z, "ZoneBC_t", 1, "BC_t", ibc, "end") == CG_OK) {
            if (cg_famname_read(famname) == CG_OK)
                entry.family = famname;
        }

        entries.push_back(entry);
    }

    for (int ibc=1; ibc<=nBC; ++ibc)
    {
        const BCEntry& entry = entries[static_cast<size_t>(ibc - 1)];
        std::string label = entry.family.empty() ? entry.name : entry.family;

        std::cout << "      • " << label;
        if (entry.type != CGNS_ENUMV(FamilySpecified))
            std::cout << "  (" << bc_type_string(entry.type) << ")";

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

/*------------------------------------------------------------*/
int main(int argc,char** argv)
{
    if (argc!=2){ std::cerr << "Usage: bc_dump mesh.cgns\n"; return 1; }

    const std::string path = argv[1];
    fs::Logger log("bc_dump.log");
    log.info("Opening "+path);

    try{
        fs::Mesh mesh; mesh.open(path,false);

        /* --- map native zone names to Zone<N> strings --- */
        //std::unordered_map<std::string,std::string> zname2nice;
        g_zname2nice.clear();
        for (const auto& z : mesh.zones()) {
            g_zname2nice[z.name] = "Z" + std::to_string(z.idx);
        }

        g_zoneNk.clear();
        for (const auto& z : mesh.zones())
            g_zoneNk.push_back(z.nk());

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
