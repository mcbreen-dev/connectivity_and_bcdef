/*─────────────────────────────────────────────────────────────
  File: src/connectivity_writer.cpp

  CGNS GridConnectivity1to1_t writer.

  This module is responsible for writing a single structured 1-to-1
  connectivity patch into a CGNS file, as a GridConnectivity1to1_t node.

  Inputs:
    - Mesh& mesh         : open CGNS file + zone metadata
    - ConnPatch cp       : resolved interface patch produced by the
                           connectivity pipeline (connectivity_verify.cpp)

  Two encoding paths are supported:
    1) Explicit-range path (preferred; cp.useExplicitRanges == true)
       - recvRange and donorRange are already computed and stored in cp
       - transform[3] is provided directly in cp.transform

    2) Legacy offset/size/orientation path (cp.useExplicitRanges == false)
       - patch described via face-local offsets/sizes and ori bits
       - this file derives receiver/donor PointRange and transform[3]

  Safety/robustness features:
    - 2D zones (nk==1) are handled by collapsing K indices and setting
      transform[2] = 0
    - existing connections are indexed and checked to avoid duplicating
      an identical GridConnectivity1to1_t entry on re-runs
    - name collision handling: if a candidate name already exists but
      points to a different range, a suffix is added (_1, _2, ...)

  Notes on CGNS conventions:
    - Structured 1-to-1 connectivity is stored per-zone under ZoneGridConnectivity_t
    - The Transform vector is a signed axis permutation describing donor
      index directions relative to receiver index directions
    - Interface ranges are expressed in full IndexDimension (I,J,K), even
      when a patch lies on a constant-index face; one component is constant
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

/*=====================================================================
  Face-axis helpers

  The connectivity pipeline uses FaceDir to represent structured block
  faces (iMin/iMax/jMin/jMax/kMin/kMax). Several helper routines in this
  file convert FaceDir into:
    - the primary axis (I/J/K) of the face normal
    - whether the face is the "plus" side (max) or "minus" side (min)
    - the two varying (in-plane) axes for that face

  These helpers are used when reconstructing:
    - the CGNS Transform vector (tvec[3])
    - point ranges for receiver/donor patches
=====================================================================*/
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

/*------------------------------------------------------------------
  face_variable_axes

  For a given face, return the two axes that vary on that face plane.

  Output:
    out[0], out[1] are axis ids (0=I, 1=J, 2=K) that span the face.
------------------------------------------------------------------*/
static void face_variable_axes(FaceDir f,int out[2])
{
    switch (f){
        case FaceDir::IMIN: case FaceDir::IMAX: out[0]=1; out[1]=2; break;
        case FaceDir::JMIN: case FaceDir::JMAX: out[0]=0; out[1]=2; break;
        case FaceDir::KMIN: case FaceDir::KMAX: out[0]=0; out[1]=1; break;
    }
}

/*=====================================================================
  axis_map

  Construct one component of the CGNS Transform vector.

  Transform is stored as signed axes:
      ±1 => donor I
      ±2 => donor J
      ±3 => donor K
      0  => collapsed axis (used for 2D)

  Parameters:
    recv, donor : receiver and donor faces (FaceDir)
    rAx         : receiver axis being mapped (0=I, 1=J, 2=K)
    swapUV      : whether receiver U/V were swapped
    flipU/V     : whether U/V axis direction was flipped

  Behavior:
    - If rAx is the collapsed axis of the receiver face (normal direction),
      map it to the donor's normal axis with sign based on min/max pairing.
    - Otherwise map receiver face tangential axes according to swap/flip.
=====================================================================*/
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


/*
CGNS-SIDS (Multizone Interface Connectivity):
Interface ranges (PointRange / PointRangeDonor) are expressed in the grid IndexDimension,
even when describing a lower-dimensional interface like a face. For structured grids this
means one index component is constant (redundant) for an abutting face patch.
This code always constructs 3-axis ranges (I,J,K) and collapses the unused axis for 2-D.

CGNS User Guide. CGNS/SIDS - Standard Interface Data Structures/8:Multizone Interface Connectivity
https://cgns.org/standard/SIDS/multizone.html
*/

/*=====================================================================
  face_range

  Build a CGNS-style PointRange array for a face patch.

  Inputs:
    z          : zone dimensions
    d          : which face of the zone
    offU/offV  : offsets along the two varying face axes
    nu/nv      : extents (in vertices) along the two varying face axes

  Output:
    rng[6] = { imin, jmin, kmin, imax, jmax, kmax } (inclusive)

  Notes:
    - Offsets are applied in vertex coordinates (not cells).
    - Callers often normalise the range after construction.
=====================================================================*/
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

/*=====================================================================
  uv_dims_local

  Return face-local vertex extents (U,V) for a given zone face.

  For a face plane:
    - i faces vary over (J,K)
    - j faces vary over (I,K)
    - k faces vary over (I,J)

  Used by the legacy (non-explicit) path to derive default full-face
  ranges and clamp patch sizes.
=====================================================================*/
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

/*=====================================================================
  Duplicate-avoidance cache

  CGNS permits multiple GridConnectivity1to1_t entries per zone. When this
  program is run repeatedly (or when overwrite is not enabled), an
  identical interface might already exist.

  To avoid duplicate entries:
    - the first time we write for a receiver zone, we read all existing
      1to1 entries and store their:
        * name
        * receiver PointRange
        * donor PointRangeDonor
    - when choosing a name, we check if a candidate name exists:
        * if name exists AND ranges match exactly => skip writing
        * if name exists but ranges differ => bump suffix and retry
=====================================================================*/
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

/*--------------------------------------------------------------
  normalise_range

  Ensure begin <= end for each axis component in a 3D PointRange.
----------------------------------------------------------------*/
static void normalise_range(cgsize_t rng[6])
{
    for (int ax=0; ax<3; ++ax)
        if (rng[ax] > rng[ax+3]) std::swap(rng[ax],rng[ax+3]);
}


/*=====================================================================
  sync_collapsed_axes

  For a face patch, one receiver axis is constant (collapsed). In CGNS,
  the transform may map that collapsed receiver axis to a donor axis.

  This helper mirrors receiver constant-axis values into donor range so
  both ranges are consistent (same constant plane index) when needed.

  Special cases:
    - if tvec[rAx] == 0, the axis is considered "collapsed because 2D"
      and no synchronization is performed for that axis.
=====================================================================*/
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


/*
CGNS-SIDS note:
A structured 1-to-1 abutting interface patch is stored twice in CGNS: one GridConnectivity1to1_t
entry under each adjacent Zone_t. (Zone-by-zone hierarchy requirement.)
This tool writes only one direction per detected pair to avoid producing duplicate entries
when re-running; the complementary entry may already exist or be produced by a separate pass.

CGNS User Guide. CGNS/SIDS - Standard Interface Data Structures/8:Multizone Interface Connectivity
https://cgns.org/standard/SIDS/multizone.html
*/

/*=====================================================================
  write_1to1

  Write one GridConnectivity1to1_t entry into the receiver zone.

  Important CGNS modeling note:
    CGNS stores a structured 1-to-1 interface "twice" (once under each
    adjacent zone). This tool intentionally writes only one direction per
    detected pair at a time (see should_write() in connectivity_utils.cpp)
    to avoid duplicate reciprocal entries in repeated runs.

  The receiver/donor zone IDs in ConnPatch are 1-based.
=====================================================================*/
void write_1to1(Mesh& mesh,const ConnPatch& cp)
{
    const Zone& recvZ = mesh.zones()[cp.recvZone-1];
    const Zone& donorZ = mesh.zones()[cp.donorZone-1];

    cgsize_t rRange[6]{};
    cgsize_t dRange[6]{};
    int tvec[3]{0, 0, 0};

    /*=========================================================
      Path A: explicit ranges

      resolve_matches_*() already computed:
        - receiver PointRange begin/end
        - donor PointRangeDonor begin/end
        - transform vector (signed axes)

      This is the preferred path in the current pipeline.
    =========================================================*/
    if (cp.useExplicitRanges) {
        tvec[0] = cp.transform[0];
        tvec[1] = cp.transform[1];
        tvec[2] = cp.transform[2];

        /*-----------------------------------------------------
          2D receiver zones are represented by nk==1.
          CGNS transform collapses K by setting transform[2]=0 and
          ranges become 2D-style (imin jmin imax jmax 1 1).
        -----------------------------------------------------*/
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

            /* collapsed axis in 2D */
            tvec[2] = 0;

            /* ensure remaining transform entries are valid non-K axes */
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
    /*=========================================================
      Path B: legacy offset/size/orientation encoding

      This path reconstructs rRange/dRange and tvec from:
        - receiver/donor face ids (FaceDir)
        - patch offsets and extents (offU/offV/sizeU/sizeV)
        - orientation bitfield cp.ori:
            bit0 flipU
            bit1 flipV
            bit2 swapUV
    =========================================================*/

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

    /*--------------------------------------------------------------
      Additional sign correction when receiver and donor are faces on
      the same axis plane (I↔I, J↔J, K↔K). This enforces CGNS sign
      conventions for min/min, min/max, max/min, max/max pairings.
    --------------------------------------------------------------*/
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

    /*=========================================================
      Existing-entry indexing for this receiver zone

      Read existing 1-to-1 entries once per receiver zone, store
      name -> ranges, and use this info to skip duplicates.
    =========================================================*/
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

                /*-----------------------------------------------------
                  Store ranges in a consistent 6-element layout.

                  For 2D zones, CGNS returns 4 values for ranges but this
                  writer always reasons in 6 values (with K collapsed to 1).
                -----------------------------------------------------*/
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

    /*=========================================================
      Name selection + duplicate guard

      Base name scheme matches bc_dump parsing conventions:
        Z<recv>_<recvFace>_to_Z<donor>_<donorFace>[_N]

      If a candidate name already exists:
        - if ranges match => skip (already written)
        - else            => increment suffix and retry
    =========================================================*/
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

            /*------------------------------------------------------
              If name exists and ranges are identical, no-op.
              This is the main "skip duplicates on re-run" behavior.
            ------------------------------------------------------*/
            if (ranges_equal(entry.recv, {rRange[0], rRange[1], rRange[2], rRange[3], rRange[4], rRange[5]}) &&
                ranges_equal(entry.donor, {dRange[0], dRange[1], dRange[2], dRange[3], dRange[4], dRange[5]})) {
                return;
            }

            /* name collision with different geometry => pick a new name */
            ++count;
            continue;
        }
        nameStr = candidate;
        ++count;
        break;
    }
    std::snprintf(name, sizeof(name), "%s", nameStr.c_str());

    /*=========================================================
      CGNS write call

      cg_1to1_write writes a GridConnectivity1to1_t node under:
        /Base/Zone[recvZone]/ZoneGridConnectivity_t

      Arguments:
        - receiver zone id
        - connection name
        - donor zone name (string)
        - receiver PointRange (IJK begin/end)
        - donor PointRangeDonor
        - transform vector (signed axes)
    =========================================================*/
    int connID=0;
    int ier = cg_1to1_write(mesh.file_id(),mesh.base_id(),
                            cp.recvZone,name,
                            donorZ.name.c_str(),
                            rRange,dRange,tvec,&connID);
    if (ier){
        cg_error_print();
        throw std::runtime_error("cg_1to1_write failed for "+std::string(name));
    }

    /*--------------------------------------------------------------
      Update existing-entry cache so additional writes in this run
      also skip duplicates and avoid name reuse.
    --------------------------------------------------------------*/
    ExistingConn new_entry;
    new_entry.name = nameStr;
    new_entry.recv = {rRange[0], rRange[1], rRange[2], rRange[3], rRange[4], rRange[5]};
    new_entry.donor = {dRange[0], dRange[1], dRange[2], dRange[3], dRange[4], dRange[5]};
    existing.by_name[new_entry.name] = existing.entries.size();
    existing.entries.push_back(new_entry);
}

} // namespace fs
