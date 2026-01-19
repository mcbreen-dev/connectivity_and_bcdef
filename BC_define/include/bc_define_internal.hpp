/*
  File: include/bc_define_internal.hpp

  Internal declarations for the bc_define executable and Plot3D text outputs.

  This header is included by:
    - src/bc_define_helpers.cpp
    - src/bc_define_bcspec.cpp
    - src/bc_define_plot3d.cpp
    - src/bc_define_cgns.cpp
    - src/bc_define_main.cpp

  It centralizes:
    - types used to key BC specifications (ZoneFaceKey) and store parsed specs (BCSpec)
    - Plot3D BC output row representation (Plot3DBcEntry)
    - parsing helpers for the bcdef.input file
    - helper functions that convert vertex-range boundary patches into face-center ranges
    - CGNS helpers used when writing BCs into a CGNS file (family creation, face detection)
    - APIs used by main() to write boundary conditions / plot3d text outputs / IO benchmark

  Notes on indexing conventions used throughout bc_define:
    - "zone" is 1-based (matches CGNS zone numbering and Plot3DZone::idx assignment).
    - FaceDir is the logical face selector (i-/i+/j-/j+/k-/k+).
    - BoundaryPatch::vtxBegin / vtxEnd are vertex indices (1-based), coming from connectivity.
    - BC writing for CGNS uses face-center ranges (cell-centered on boundary faces):
        * for a vertex range [lo..hi], the corresponding face-center range is [lo..hi-1]
      with special handling for 2D collapsed dimensions.

  Size / type notes:
    - ZoneFaceKey stores zone as int (typically 32-bit) and face as FaceDir (uint8_t).
    - Many CGNS indices use cgsize_t; this header frequently converts to long long
      or cgsize_t arrays depending on the call site.
*/
#pragma once

#include "common.hpp"
#include "mesh_io.hpp"
#include "plot3d_io.hpp"
#include "connectivity.hpp"

#include <algorithm>
#include <array>
#include <string>
#include <unordered_map>
#include <vector>

#include <cgnslib.h>

#ifndef CGNS_MAX_NAME_LENGTH
#define CGNS_MAX_NAME_LENGTH 32
#endif

namespace bcdef {
class Logger;
}

namespace bcdef::boundary {

/*
  ZoneFaceKey

  Key type used to associate a BC specification with a specific zone face.

  Used in:
    - parse_bcdef(): returns unordered_map<ZoneFaceKey, BCSpec, ZoneFaceKeyHash>
    - write_boundary_conditions(): finds a spec for each detected BoundaryPatch
    - build_plot3d_bcs(): same mapping for Plot3D outputs

  Fields:
    zone:
      - 1-based zone index (int)
      - for CGNS: equals Zone::idx
      - for Plot3D: equals Plot3DZone::idx
    face:
      - FaceDir enum (stored as uint8_t in common.hpp)

  Equality:
    - both zone and face must match.
*/
struct ZoneFaceKey {
    int zone = 0;
    bcdef::FaceDir face = bcdef::FaceDir::IMIN;

    bool operator==(const ZoneFaceKey& other) const
    {
        return zone == other.zone && face == other.face;
    }
};

/*
  ZoneFaceKeyHash

  Hash functor for ZoneFaceKey for unordered_map usage.

  Implementation:
    - shifts zone by 3 bits and XORs with face value
    - face has 6 values (0..5), fits in the lower 3 bits.
*/
struct ZoneFaceKeyHash {
    size_t operator()(const ZoneFaceKey& key) const
    {
        return (static_cast<size_t>(key.zone) << 3) ^
               static_cast<size_t>(key.face);
    }
};

/*
  BCSpec

  Parsed BC assignment for a zone face.

  Produced by:
    - parse_bcdef(path)

  Used by:
    - write_boundary_conditions() (CGNS writing)
    - build_plot3d_bcs() (Plot3D text output)

  Fields:
    family:
      - the Family_t name (string) used when writing CGNS BCs as FamilySpecified
      - for Plot3D output: written as the "family" column in "bcs" file

    type:
      - CGNS_ENUMT(BCType_t)
      - default: BCTypeUserDefined if no mapping is found
      - parse_bc_type() maps common tokens (wall, farfield, inflow, etc.) to BCType_t
*/
struct BCSpec {
    std::string family;
    CGNS_ENUMT(BCType_t) type = CGNS_ENUMV(BCTypeUserDefined);
};

/*
  Plot3DBcEntry

  One row in the Plot3D "bcs" output file produced by write_plot3d_bcs().

  Produced by:
    - build_plot3d_bcs(), then write_plot3d_bcs()

  Fields:
    zone:
      - 1-based zone index (int)

    face:
      - FaceDir selector for which outer face this range lies on

    family:
      - BC "family" label written into the output file

    range:
      - 6 values: [i1,j1,k1,i2,j2,k2] as long long
      - this is a face-center range (cell-space on the face), not vertex-space
      - produced by face_center_range(zone, face, vtxBegin, vtxEnd)
*/
struct Plot3DBcEntry {
    int zone = 0;
    bcdef::FaceDir face = bcdef::FaceDir::IMIN;
    std::string family;
    std::array<long long, 6> range{};
};

/*------------------------------------------------------------------------------
  Text parsing helpers (src/bc_define_helpers.cpp, src/bc_define_bcspec.cpp)
------------------------------------------------------------------------------*/

/*
  trim(s)

  Removes leading and trailing ASCII whitespace using std::isspace.
  Used when parsing bcdef.input lines and comma-separated zone lists.
*/
std::string trim(const std::string& s);

/*
  normalize_token(s)

  Lowercases and keeps only alphanumeric characters.
  Used to compare BC tokens in a case/format-insensitive way.
*/
std::string normalize_token(const std::string& s);

/*
  parse_zone_list(token)

  Parses zone specifications like:
    "1"          -> {1}
    "1,2,5"      -> {1,2,5}
    "3-7"        -> {3,4,5,6,7}
    "1-3,7,9-10" -> {1,2,3,7,9,10}

  Returns:
    - sorted unique vector<int>

  Used by:
    - parse_bcdef()
*/
std::vector<int> parse_zone_list(const std::string& token);

/*
  parse_bc_type(token)

  Converts a user token (e.g. "wall", "farfield", "symmetry") into a CGNS BCType_t.
  Implementation is in src/bc_define_bcspec.cpp.

  Used by:
    - parse_bcdef()
*/
CGNS_ENUMT(BCType_t) parse_bc_type(const std::string& token);

/*
  is_plot3d_mesh(path)

  Returns true if filesystem extension is ".x" (case-insensitive).
  Used in bc_define_main.cpp to pick Plot3D vs CGNS workflow.
*/
bool is_plot3d_mesh(const std::string& path);

/*------------------------------------------------------------------------------
  Range conversion helpers (vertex-space -> face-center-space)
------------------------------------------------------------------------------*/

/*
  face_center_range(zone, face, vtxBegin, vtxEnd)  (CGNS Zone overload)

  Input:
    - zone: bcdef::Zone containing Ni/Nj/Nk vertex dimensions
    - face: which boundary face the patch lies on
    - vtxBegin, vtxEnd:
        vertex indices (1-based) marking a patch on the face

  Output:
    - bcdef::PointRange in face-center (cell-face) indexing:
        For non-collapsed dimensions: [lo .. hi-1]
        For collapsed dimensions (dim<=1): [1..1]

  Used in:
    - write_boundary_conditions() (to compute cgsize_t ranges for cg_boco_write)
*/
bcdef::PointRange face_center_range(const bcdef::Zone& zone,
                                 bcdef::FaceDir face,
                                 const bcdef::IJK& vtxBegin,
                                 const bcdef::IJK& vtxEnd);
/*
  face_center_range(zone, face, vtxBegin, vtxEnd)  (Plot3DZone overload)

  Same logic as the Zone overload, but uses Plot3DZone::ni/nj/nk.

  Used in:
    - build_plot3d_bcs()
*/
bcdef::PointRange face_center_range(const bcdef::Plot3DZone& zone,
                                 bcdef::FaceDir face,
                                 const bcdef::IJK& vtxBegin,
                                 const bcdef::IJK& vtxEnd);

/*------------------------------------------------------------------------------
  CGNS writing helpers (src/bc_define_cgns.cpp)
------------------------------------------------------------------------------*/

/*
  face_grid_location(face, cell_dim)

  Returns the CGNS GridLocation_t appropriate for a boundary patch on a face.

  Behavior (src/bc_define_cgns.cpp):
    - cell_dim >= 3:
        IMIN/IMAX -> IFaceCenter
        JMIN/JMAX -> JFaceCenter
        KMIN/KMAX -> KFaceCenter
    - cell_dim == 2:
        -> EdgeCenter
    - otherwise: throws

  Used by:
    - write_boundary_conditions() to set grid location for each new BC.
*/
CGNS_ENUMT(GridLocation_t) face_grid_location(bcdef::FaceDir face, int cell_dim);

/*
  Template: ranges_equal(a,b)

  Exact equality comparison for fixed-size 6-element ranges.
  Used for:
    - skipping rewrites when the same range already exists
    - overlap handling in both CGNS and Plot3D workflows
*/
template <typename T>
inline bool ranges_equal(const std::array<T, 6>& a, const std::array<T, 6>& b)
{
    return a == b;
}

/*
  Template: range_overlap(a,b,out)

  Computes overlap between two axis-aligned inclusive index ranges:
    a = [a0,a1,a2,a3,a4,a5]
    b = [b0,b1,b2,b3,b4,b5]
  where dims 0..2 are mins and dims 3..5 are maxes.

  Output:
    - out receives the overlap range if overlap exists

  Returns:
    - true if overlap exists (including touching boundary, since lo>hi is the fail case)

  Used in:
    - write_boundary_conditions(): detects overlaps against existing BCs
    - build_plot3d_bcs(): detects overlaps across emitted entries
*/
template <typename T>
inline bool range_overlap(const std::array<T, 6>& a, const std::array<T, 6>& b,
                          std::array<T, 6>& out)
{
    for (int dim = 0; dim < 3; ++dim) {
        T lo = std::max(a[dim], b[dim]);
        T hi = std::min(a[dim + 3], b[dim + 3]);
        if (lo > hi)
            return false;
        out[dim] = lo;
        out[dim + 3] = hi;
    }
    return true;
}

/*
  Template: range_to_string(r)

  Debug/exception formatting helper for 6-element index ranges.
  Used in overlap error messages in bc_define_cgns.cpp and bc_define_plot3d.cpp.
*/
template <typename T>
inline std::string range_to_string(const std::array<T, 6>& r)
{
    return "[" + std::to_string(r[0]) + "," + std::to_string(r[1]) + "," +
           std::to_string(r[2]) + "]-[" + std::to_string(r[3]) + "," +
           std::to_string(r[4]) + "," + std::to_string(r[5]) + "]";
}

/*
  bc_patch_name(family, face, index)

  Builds a BC_t node name:
    "<family>_<i-/i+/j-/j+/k-/k+>_<index>"

  Behavior:
    - truncates to CGNS_MAX_NAME_LENGTH (32) to satisfy CGNS name constraints.

  Used in:
    - write_boundary_conditions() when creating new BC_t nodes.
*/
std::string bc_patch_name(const std::string& family, bcdef::FaceDir face, int index);

/*
  range_face(zone, r, out)

  Determines which outer face a face-center range lies on.

  Input:
    - zone: to know ni/nj/nk and the face-center outer indices (1 or ni-1, etc.)
    - r: 6-element face-center range (cgsize_t)
  Output:
    - out: FaceDir if r lies on an outer face
  Returns:
    - true if r corresponds to IMIN/IMAX/JMIN/JMAX or (for 3D) KMIN/KMAX

  Used in:
    - write_boundary_conditions() to compare candidate patches against existing BCs
      on the same face.
*/
bool range_face(const bcdef::Zone& zone,
                const std::array<cgsize_t, 6>& r,
                bcdef::FaceDir& out);

/*
  ensure_family(fn, B, family, bc_type, cache)

  Ensures the CGNS Base_t contains a Family_t named `family` and that the
  Family_t has a FamilyBC_t child.

  Inputs:
    fn, B:
      - CGNS file/base ids
    family:
      - Family_t name to find or create
    bc_type:
      - BCType_t to assign when creating FamilyBC if missing
    cache:
      - maps family name -> family index (int) to avoid re-scanning cg_nfamilies

  Returns:
    - Family index (1-based in CGNS) as int

  Used in:
    - write_boundary_conditions() before writing any BCs that reference a family.
*/
int ensure_family(int fn, int B,
                  const std::string& family,
                  CGNS_ENUMT(BCType_t) bc_type,
                  std::unordered_map<std::string, int>& cache);

/*------------------------------------------------------------------------------
  High-level operations (called from bc_define_main.cpp)
------------------------------------------------------------------------------*/

/*
  parse_bcdef(path)

  Reads bcdef.input and returns a map of per-zone-face BC specifications.
  Implementation: src/bc_define_bcspec.cpp

  Input line format:
    <bc_token> <zone_list> <face_selector>
  Example:
    wall 1-4 i-
    farfield 5,7,9 j+

  Output map key:
    ZoneFaceKey{zone, FaceDir}

  Duplicate handling:
    - identical repeat is allowed
    - conflicting definitions for the same zone/face throw
*/
std::unordered_map<ZoneFaceKey, BCSpec, ZoneFaceKeyHash>
parse_bcdef(const std::string& path);

/*
  write_boundary_conditions(mesh, patches, specs, autowall, autofarfield)

  Writes BC_t nodes into the CGNS file open in `mesh`.

  Inputs:
    mesh:
      - open CGNS mesh (modify mode in bc_define_main.cpp unless bench)
    patches:
      - boundary patches produced by ConnectivityDetector::run(...)
      - each patch is a face patch on a zone (possibly partial)
    specs:
      - user bcdef specs keyed by (zone, face)
    autowall/autofarfield:
      - if a boundary patch face has no user spec:
          * autowall writes BCWall family "wall"
          * autofarfield writes BCFarfield family "farfield"
        mutually exclusive in main()

  Overlap behavior:
    - checks against existing BCs in the file:
        * if user_spec overlaps and is not identical to same family/type -> throws
        * if auto_spec overlaps any existing BC -> skips
    - also skips writing if the same BC name already exists with the same range

  Implementation: src/bc_define_cgns.cpp
*/
void write_boundary_conditions(
    bcdef::Mesh& mesh,
    const std::vector<bcdef::BoundaryPatch>& patches,
    const std::unordered_map<ZoneFaceKey, BCSpec, ZoneFaceKeyHash>& specs,
    bool autowall,
    bool autofarfield);

/*
  run_io_benchmark(mesh, log, iters)

  Benchmarks CGNS coordinate reads:
    - volume reads (full CoordinateX/Y/Z)
    - face reads (CoordinateX/Y/Z on each face plane)

  Prints per-zone and overall stats to stdout and logs start/finish.
  Used only in bc_define_main.cpp when --bench-io is provided.
*/
void run_io_benchmark(bcdef::Mesh& mesh, bcdef::Logger& log, int iters);


/*
  write_plot3d_1to1s(path, zones, conns)

  Writes a text file containing all ConnPatch explicit ranges and transforms:
    recv_zone donor_zone
    recv range (i1 j1 k1 i2 j2 k2)
    donor range (i1 j1 k1 i2 j2 k2)
    transform (t1 t2 t3)

  Called in bc_define_main.cpp (Plot3D branch):
    - writes file "1to1s" in current working directory.
*/
void write_plot3d_1to1s(const std::string& path,
                        const std::vector<bcdef::Plot3DZone>& zones,
                        const std::vector<bcdef::ConnPatch>& conns);

/*
  build_plot3d_bcs(zones, patches, specs, autowall, autofarfield)

  Builds Plot3DBcEntry list for boundary patches in Plot3D workflow.

  Behavior:
    - for each BoundaryPatch:
        * looks up ZoneFaceKey in specs
        * if missing and auto flags are set, uses "wall" or "farfield"
        * converts patch vertex range -> face-center range
        * prevents overlaps within a zone/face:
            - user spec overlap throws unless identical (same family and exact range)
            - auto spec overlap is skipped
    - returns a flat vector of Plot3DBcEntry in encounter order

  Called in:
    - bc_define_main.cpp (Plot3D branch)
*/
std::vector<Plot3DBcEntry> build_plot3d_bcs(
    const std::vector<bcdef::Plot3DZone>& zones,
    const std::vector<bcdef::BoundaryPatch>& patches,
    const std::unordered_map<ZoneFaceKey, BCSpec, ZoneFaceKeyHash>& specs,
    bool autowall,
    bool autofarfield);


/*
  write_plot3d_bcs(path, zones, bcs)

  Writes Plot3DBcEntry rows to a text file with header comments.
  Called in bc_define_main.cpp (Plot3D branch):
    - writes file "bcs" in current working directory.
*/
void write_plot3d_bcs(const std::string& path,
                      const std::vector<bcdef::Plot3DZone>& zones,
                      const std::vector<Plot3DBcEntry>& bcs);

} // namespace bcdef::boundary
