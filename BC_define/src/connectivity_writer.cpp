/*─────────────────────────────────────────────────────────────
  Connectivity-writer – writes one CGNS GridConnectivity1to1_t
─────────────────────────────────────────────────────────────*/
#include "connectivity_writer.hpp"
#include "common.hpp"
#include "logger.hpp"
#include "mesh_io.hpp"

#include <array>
#include <unordered_map>
#include <string>
#include <cstring>
#include <iostream>

#ifndef CGNS_MAX_NAME_LENGTH
#define CGNS_MAX_NAME_LENGTH 32
#endif

namespace fs {

/* same 3-axis table we use for quick 2-D swaps/flips */
const int map[8][3] = {
    { 1,  2,  3}, {-1,  2,  3}, { 1, -2,  3}, {-1, -2,  3},
    { 2,  1, -3}, {-2,  1, -3}, { 2, -1, -3}, {-2, -1, -3}
};

/*--------------------------------------------------------------
   helpers (axis index, plus-side?, etc.)
----------------------------------------------------------------*/
static constexpr int face_axis(FaceDir f)
{
    switch (f) {
        case FaceDir::IMIN: case FaceDir::IMAX: return 0;
        case FaceDir::JMIN: case FaceDir::JMAX: return 1;
        case FaceDir::KMIN: case FaceDir::KMAX: return 2;
    }
    return 0;
}
static constexpr bool is_plus_side(FaceDir f)
{
    return f==FaceDir::IMAX || f==FaceDir::JMAX || f==FaceDir::KMAX;
}

/*------------------------------------------------------------------*/
static void face_variable_axes(FaceDir f,int out[2])
{
    switch (f){
        case FaceDir::IMIN: case FaceDir::IMAX: out[0]=1; out[1]=2; break;
        case FaceDir::JMIN: case FaceDir::JMAX: out[0]=0; out[1]=2; break;
        case FaceDir::KMIN: case FaceDir::KMAX: out[0]=0; out[1]=1; break;
    }
}

/*------------------------------------------------------------------*/
static int axis_map(FaceDir recv,FaceDir donor,int rAx,
                    bool swapUV,bool flipU,bool flipV)
{
    int recvVar[2], donorVar[2];
    face_variable_axes(recv, recvVar);
    face_variable_axes(donor,donorVar);

    if (rAx == face_axis(recv)){                 /* collapsed axis */
        int sign = is_plus_side(recv)==is_plus_side(donor) ? -1:+1;
        return sign * (face_axis(donor)+1);      /* ±(1|2|3) */
    }

    int idx = (rAx == recvVar[swapUV?1:0]) ? 0 : 1;   /* which param (U/V) */
    int dAx = donorVar[idx];
    int sign = (idx==0 ? (flipU?-1:+1) : (flipV?-1:+1));
    return sign * (dAx+1);
}

/*------------------------------------------------------------------*/
static void face_range(const Zone& z,FaceDir d,
                       long long offU,long long offV,
                       long long nu,long long nv,
                       cgsize_t rng[6])     /* imin jmin kmin imax jmax kmax */
{
    long long Ni=z.ni(), Nj=z.nj(), Nk=z.nk();
    long long imin=1,imax=Ni, jmin=1,jmax=Nj,kmin=1,kmax=Nk;

    switch (d){
        case FaceDir::IMIN: imin=imax=1;
            jmin+=offU; jmax=offU+nu;   kmin+=offV; kmax=offV+nv; break;
        case FaceDir::IMAX: imin=imax=Ni;
            jmin+=offU; jmax=offU+nu;   kmin+=offV; kmax=offV+nv; break;

        case FaceDir::JMIN: jmin=jmax=1;
            imin+=offU; imax=offU+nu;   kmin+=offV; kmax=offV+nv; break;
        case FaceDir::JMAX: jmin=jmax=Nj;
            imin+=offU; imax=offU+nu;   kmin+=offV; kmax=offV+nv; break;

        case FaceDir::KMIN: kmin=kmax=1;
            imin+=offU; imax=offU+nu;   jmin+=offV; jmax=offV+nv; break;
        case FaceDir::KMAX: kmin=kmax=Nk;
            imin+=offU; imax=offU+nu;   jmin+=offV; jmax=offV+nv; break;
    }
    rng[0]=imin; rng[1]=jmin; rng[2]=kmin;
    rng[3]=imax; rng[4]=jmax; rng[5]=kmax;
}

/*------------------------------------------------------------------*/
static inline std::pair<long long,long long>
uv_dims_local(const Zone& z,FaceDir d)
{
    switch (d){
        case FaceDir::IMIN: case FaceDir::IMAX: return {z.nj(), z.nk()};
        case FaceDir::JMIN: case FaceDir::JMAX: return {z.ni(), z.nk()};
        case FaceDir::KMIN: case FaceDir::KMAX: return {z.ni(), z.nj()};
    }
    return {0,0};
}

/*------------------------------------------------------------------*/
namespace {
    struct ExistingConn {
        std::string name;
        std::array<cgsize_t, 6> recv{};
        std::array<cgsize_t, 6> donor{};
    };
    struct ZoneConnIndex {
        std::unordered_map<std::string, size_t> by_name;
        std::vector<ExistingConn> entries;
        bool loaded = false;
    };

    std::unordered_map<std::string, int> g_name_counts;
    std::unordered_map<int, ZoneConnIndex> g_existing_conns;
}

/*------------------------------------------------------------------*/
static void normalise_range(cgsize_t rng[6])
{
    for (int ax=0; ax<3; ++ax)
        if (rng[ax] > rng[ax+3]) std::swap(rng[ax],rng[ax+3]);
}

/*------------------------------------------------------------------*/
static void sync_collapsed_axes(const cgsize_t r[6],
                                cgsize_t       d[6],
                                const int      tvec[3])
{
    for (int rAx=0; rAx<3; ++rAx){
        if (r[rAx]!=r[rAx+3]) continue;         /* varies */
        if (tvec[rAx] == 0) continue;           /* 2-D collapsed axis */
        int dAx = std::abs(tvec[rAx])-1;
        d[dAx]   = r[rAx];
        d[dAx+3] = r[rAx+3];
    }
}

/*====================================================================*/
void write_1to1(Mesh& mesh,const ConnPatch& cp)
{
    const Zone& recvZ = mesh.zones()[cp.recvZone-1];
    const Zone& donorZ = mesh.zones()[cp.donorZone-1];

    cgsize_t rRange[6]{};
    cgsize_t dRange[6]{};
    int tvec[3]{0, 0, 0};

    if (cp.useExplicitRanges) {
        tvec[0] = cp.transform[0];
        tvec[1] = cp.transform[1];
        tvec[2] = cp.transform[2];

        if (recvZ.nk() == 1) {
            rRange[0] = static_cast<cgsize_t>(cp.recvRange.begin[0]);
            rRange[1] = static_cast<cgsize_t>(cp.recvRange.begin[1]);
            rRange[2] = static_cast<cgsize_t>(cp.recvRange.end[0]);
            rRange[3] = static_cast<cgsize_t>(cp.recvRange.end[1]);
            rRange[4] = 1;
            rRange[5] = 1;

            dRange[0] = static_cast<cgsize_t>(cp.donorRange.begin[0]);
            dRange[1] = static_cast<cgsize_t>(cp.donorRange.begin[1]);
            dRange[2] = static_cast<cgsize_t>(cp.donorRange.end[0]);
            dRange[3] = static_cast<cgsize_t>(cp.donorRange.end[1]);
            dRange[4] = 1;
            dRange[5] = 1;
            tvec[2] = 0;
            for (int k = 0; k < 2; ++k) {
                if (tvec[k] == 0 || std::abs(tvec[k]) == 3)
                    tvec[k] = (k == 0 ? 1 : 2);
            }
        } else {
            rRange[0] = static_cast<cgsize_t>(cp.recvRange.begin[0]);
            rRange[1] = static_cast<cgsize_t>(cp.recvRange.begin[1]);
            rRange[2] = static_cast<cgsize_t>(cp.recvRange.begin[2]);
            rRange[3] = static_cast<cgsize_t>(cp.recvRange.end[0]);
            rRange[4] = static_cast<cgsize_t>(cp.recvRange.end[1]);
            rRange[5] = static_cast<cgsize_t>(cp.recvRange.end[2]);

            dRange[0] = static_cast<cgsize_t>(cp.donorRange.begin[0]);
            dRange[1] = static_cast<cgsize_t>(cp.donorRange.begin[1]);
            dRange[2] = static_cast<cgsize_t>(cp.donorRange.begin[2]);
            dRange[3] = static_cast<cgsize_t>(cp.donorRange.end[0]);
            dRange[4] = static_cast<cgsize_t>(cp.donorRange.end[1]);
            dRange[5] = static_cast<cgsize_t>(cp.donorRange.end[2]);
        }
    } else {
    /* -- receiver PointRange --------------------------------------- */
    auto [fullU,fullV] = uv_dims_local(recvZ,cp.recvDir);

    long long offU = cp.offU, offV = cp.offV;
    if (offU<0 || offU>=fullU) offU = 0;
    if (offV<0 || offV>=fullV) offV = 0;

    long long sizeU = cp.sizeU<=0 || cp.sizeU>fullU ? fullU : cp.sizeU;
    long long sizeV = cp.sizeV<=0 || cp.sizeV>fullV ? fullV : cp.sizeV;

    face_range(recvZ,cp.recvDir,offU,offV,sizeU,sizeV,rRange);
    normalise_range(rRange);

    /* -- derive final tvec from ori bits --------------------------- */
    {
        bool swapUV = cp.ori & 4;
        bool flipU  = cp.ori & 1;
        bool flipV  = cp.ori & 2;
        tvec[0] = axis_map(cp.recvDir,cp.donorDir,0,swapUV,flipU,flipV);
        tvec[1] = axis_map(cp.recvDir,cp.donorDir,1,swapUV,flipU,flipV);
        tvec[2] = axis_map(cp.recvDir,cp.donorDir,2,swapUV,flipU,flipV);
    }

    {
        int ax = face_axis(cp.recvDir);        // 0=I,1=J,2=K
        if (ax == face_axis(cp.donorDir)) {    // both faces lie on same axis plane
            bool same_side = (is_plus_side(cp.recvDir) == is_plus_side(cp.donorDir));
            int sgn = same_side ? +1 : -1;
            tvec[ax] =  sgn * std::abs(tvec[ax]);   // enforce correct sign
            if (same_side && ax == 0 && cp.ori == 1) {
                tvec[0] = -tvec[0];
                tvec[1] = -tvec[1];
                tvec[2] = -tvec[2];
            }
        }
    }

    /* 2-D safety (collapsed K) */
    if (recvZ.nk()==1) {
        tvec[2]=0;
        for (int k=0;k<2;++k)
            if (tvec[k]==0 || std::abs(tvec[k])==3)
                tvec[k] = (k==0?1:2);
    }

    /* -- donor PointRange (face-local) ------------------------------ */
    auto [donU,donV]   = uv_dims_local(donorZ,cp.donorDir);

    face_range(donorZ,cp.donorDir,0,0,donU,donV,dRange);
    normalise_range(dRange);

    /* apply patch offset (if any) */
    if (cp.offU || cp.offV){
        int dU = std::abs(tvec[0])-1;
        int dV = std::abs(tvec[1])-1;
        dRange[dU]     += cp.offU; dRange[dU+3] += cp.offU;
        dRange[dV]     += cp.offV; dRange[dV+3] += cp.offV;
    }

    /* mirror collapsed axes */
    sync_collapsed_axes(rRange,dRange,tvec);

    /* -- 2-D remap for CGNS (index_dim=2) --------------------------- */
    if (recvZ.nk() == 1) {
        cgsize_t r2[6] = {rRange[0], rRange[1], rRange[3], rRange[4], 1, 1};
        cgsize_t d2[6] = {dRange[0], dRange[1], dRange[3], dRange[4], 1, 1};
        for (int i = 0; i < 6; ++i) {
            rRange[i] = r2[i];
            dRange[i] = d2[i];
        }
    }
    }

    auto ranges_equal = [](const std::array<cgsize_t, 6>& a,
                           const std::array<cgsize_t, 6>& b) {
        return a == b;
    };

    auto& existing = g_existing_conns[cp.recvZone];
    if (!existing.loaded) {
        int nLink = 0;
        if (cg_n1to1(mesh.file_id(), mesh.base_id(), cp.recvZone, &nLink) == CG_OK) {
            for (int i = 1; i <= nLink; ++i) {
                char name[CGNS_MAX_NAME_LENGTH + 1] = {};
                char donor[CGNS_MAX_NAME_LENGTH + 1] = {};
                cgsize_t range[6] = {0};
                cgsize_t donor_range[6] = {0};
                int T[3] = {0, 0, 0};
                if (cg_1to1_read(mesh.file_id(), mesh.base_id(), cp.recvZone, i,
                                 name, donor, range, donor_range, T) != CG_OK)
                    continue;

                ExistingConn entry;
                entry.name = name;
                if (recvZ.nk() == 1) {
                    entry.recv = {range[0], range[1], range[2], range[3], 1, 1};
                    entry.donor = {donor_range[0], donor_range[1], donor_range[2], donor_range[3], 1, 1};
                } else {
                    entry.recv = {range[0], range[1], range[2], range[3], range[4], range[5]};
                    entry.donor = {donor_range[0], donor_range[1], donor_range[2],
                                   donor_range[3], donor_range[4], donor_range[5]};
                }
                existing.by_name[entry.name] = existing.entries.size();
                existing.entries.push_back(entry);
            }
        }
        existing.loaded = true;
    }

    /* -- write, guarding against duplicates ------------------------ */
    char name[64];
    std::string base = "Z" + std::to_string(cp.recvZone) + "_" +
                       facedir_to_string(cp.recvDir) + "_to_Z" +
                       std::to_string(cp.donorZone) + "_" +
                       facedir_to_string(cp.donorDir);
    int& count = g_name_counts[base];
    std::string nameStr;
    while (true) {
        std::string candidate = (count == 0) ? base : base + "_" + std::to_string(count);
        auto it = existing.by_name.find(candidate);
        if (it != existing.by_name.end()) {
            const ExistingConn& entry = existing.entries[it->second];
            if (ranges_equal(entry.recv, {rRange[0], rRange[1], rRange[2], rRange[3], rRange[4], rRange[5]}) &&
                ranges_equal(entry.donor, {dRange[0], dRange[1], dRange[2], dRange[3], dRange[4], dRange[5]})) {
                return;
            }
            ++count;
            continue;
        }
        nameStr = candidate;
        ++count;
        break;
    }
    std::snprintf(name, sizeof(name), "%s", nameStr.c_str());

    int connID=0;
    int ier = cg_1to1_write(mesh.file_id(),mesh.base_id(),
                            cp.recvZone,name,
                            donorZ.name.c_str(),
                            rRange,dRange,tvec,&connID);
    if (ier){
        cg_error_print();
        throw std::runtime_error("cg_1to1_write failed for "+std::string(name));
    }

    ExistingConn new_entry;
    new_entry.name = nameStr;
    new_entry.recv = {rRange[0], rRange[1], rRange[2], rRange[3], rRange[4], rRange[5]};
    new_entry.donor = {dRange[0], dRange[1], dRange[2], dRange[3], dRange[4], dRange[5]};
    existing.by_name[new_entry.name] = existing.entries.size();
    existing.entries.push_back(new_entry);
}

} // namespace fs
