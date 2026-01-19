/*─────────────────────────────────────────────────────────────
  File: src/bc_define_helpers.cpp
  Shared small utilities for bc_define
─────────────────────────────────────────────────────────────
  This file collects small, pure helper routines used across the
  bc_define toolchain.

  Responsibilities:
    - Whitespace trimming and token normalization for config parsing
    - Zone selector parsing (single values, ranges, comma lists)
    - Simple mesh-type detection (CGNS vs Plot3D) from filename suffix
    - Stable BC patch naming with CGNS name-length constraints
    - Conversion from vertex-index patch bounds to *cell-center* ranges
      suitable for BC definitions (FaceCenter / IFaceCenter / etc.)

  Notes on index conventions:
    - User-facing and internal zone/face selection uses 1-based
      vertex indices (I,J,K) consistent with CGNS / Plot3D conventions.
    - BC point ranges in this codebase are written on *cell centers*.
      For structured grids, that means valid indices are 1..(N-1) in
      each varying direction.
    - face_center_range() performs the vertex→cell conversion by taking
      a vertex span [lo..hi] and mapping it to cell indices [lo..hi-1].
      Degenerate dimensions (dim<=1) collapse to 1..1.

  This module intentionally avoids any CGNS I/O and has no side effects
  beyond returning computed values.
─────────────────────────────────────────────────────────────*/
#include "bc_define_internal.hpp"

#include <algorithm>
#include <cctype>
#include <filesystem>

namespace bcdef::boundary {

/*=====================================================================
  trim

  Remove leading/trailing ASCII whitespace from a string.

  Used by:
    - bc_define_bcspec.cpp while parsing bcdef lines
    - parse_zone_list() while splitting comma-separated tokens
=====================================================================*/
std::string trim(const std::string& s)
{
    size_t start = 0;
    while (start < s.size() && std::isspace(static_cast<unsigned char>(s[start])))
        ++start;
    size_t end = s.size();
    while (end > start && std::isspace(static_cast<unsigned char>(s[end - 1])))
        --end;
    return s.substr(start, end - start);
}

/*=====================================================================
  normalize_token

  Canonicalize a free-form token for robust comparison.

  Normalization rules:
    - Keep only alphanumeric characters
    - Convert to lowercase

  This makes user inputs tolerant of:
    - case differences (Wall vs wall)
    - punctuation / separators (BC-Wall, bc_wall, etc.)

  Used by:
    - bc type parsing (parse_bc_type)
    - family-key mapping to preserve stable family names
=====================================================================*/
std::string normalize_token(const std::string& s)
{
    std::string out;
    out.reserve(s.size());
    for (char c : s) {
        if (std::isalnum(static_cast<unsigned char>(c)))
            out.push_back(static_cast<char>(std::tolower(static_cast<unsigned char>(c))));
    }
    return out;
}

/*=====================================================================
  parse_zone_list

  Parse a zone selector token into a sorted unique list of zone IDs.

  Supported forms (whitespace is trimmed per component):
    - "3"         -> {3}
    - "1-4"       -> {1,2,3,4}
    - "1,3,7-9"   -> {1,3,7,8,9}

  Behavior:
    - Ranges with a>b are accepted and swapped (e.g. "4-1")
    - Output is sorted and deduplicated
    - Invalid numeric fragments will throw via std::stoi

  Used by:
    - bc_define_bcspec.cpp to expand per-line zone selections
=====================================================================*/
std::vector<int> parse_zone_list(const std::string& token)
{
    std::vector<int> zones;
    size_t start = 0;
    while (start < token.size()) {
        size_t comma = token.find(',', start);
        std::string part = trim(token.substr(
            start, comma == std::string::npos ? std::string::npos : comma - start));
        if (!part.empty()) {
            size_t dash = part.find('-');
            if (dash != std::string::npos) {
                int a = std::stoi(part.substr(0, dash));
                int b = std::stoi(part.substr(dash + 1));
                if (a > b) std::swap(a, b);
                for (int v = a; v <= b; ++v)
                    zones.push_back(v);
            } else {
                zones.push_back(std::stoi(part));
            }
        }
        if (comma == std::string::npos)
            break;
        start = comma + 1;
    }

    std::sort(zones.begin(), zones.end());
    zones.erase(std::unique(zones.begin(), zones.end()), zones.end());
    return zones;
}

/*=====================================================================
  is_plot3d_mesh

  Lightweight mesh-type sniffing based on filename extension.

  Convention used by this project:
    - ".x" indicates a Plot3D multi-block grid file

  This is only a heuristic; higher-level code may still validate the file
  by attempting to parse it.
=====================================================================*/
bool is_plot3d_mesh(const std::string& path)
{
    std::string ext = std::filesystem::path(path).extension().string();
    std::string ext_lower;
    ext_lower.reserve(ext.size());
    for (char c : ext) {
        ext_lower.push_back(static_cast<char>(std::tolower(static_cast<unsigned char>(c))));
    }
    return ext_lower == ".x";
}

/*=====================================================================
  bc_patch_name

  Construct a BC_t node name for a patch written under ZoneBC_t.

  Naming scheme:
      <family>_<face>_<index>

  Notes:
    - CGNS imposes a maximum name length (CGNS_MAX_NAME_LENGTH).
    - If truncation occurs, the name is clipped to that maximum length.
      (This preserves CGNS validity but can reduce uniqueness if users
       choose extremely long family names.)
=====================================================================*/
std::string bc_patch_name(const std::string& family, bcdef::FaceDir face, int index)
{
    std::string name = family + "_" + bcdef::facedir_to_string(face) + "_" + std::to_string(index);
    if (name.size() > CGNS_MAX_NAME_LENGTH)
        name.resize(CGNS_MAX_NAME_LENGTH);
    return name;
}

/*=====================================================================
  face_center_range (Zone overload)

  Convert a vertex-index patch range on a zone face into a PointRange
  expressed in *cell-center* indices for BC writing.

  Inputs:
    - zone     : structured CGNS Zone metadata (Ni,Nj,Nk are vertex sizes)
    - face     : which face the patch lies on
    - vtxBegin : 1-based vertex indices at one corner of the patch
    - vtxEnd   : 1-based vertex indices at the opposite corner

  Output:
    - PointRange pr in (I,J,K) index space, where varying axes are cell
      indices in 1..(N-1) and the collapsed face-normal axis is fixed:
        IMIN -> i==1
        IMAX -> i==Ni-1   (cell-center plane adjacent to max vertex)
        etc.

  Key conversion:
    - For a varying axis, a vertex span [lo..hi] maps to cell indices
      [lo .. hi-1]. This corresponds to the cells between successive
      vertices.
=====================================================================*/
bcdef::PointRange face_center_range(const bcdef::Zone& zone, bcdef::FaceDir face,
                                 const bcdef::IJK& vtxBegin, const bcdef::IJK& vtxEnd)
{
    /* Normalize corner ordering so vmin <= vmax component-wise. */
    bcdef::IJK vmin = {
        std::min(vtxBegin[0], vtxEnd[0]),
        std::min(vtxBegin[1], vtxEnd[1]),
        std::min(vtxBegin[2], vtxEnd[2])
    };
    bcdef::IJK vmax = {
        std::max(vtxBegin[0], vtxEnd[0]),
        std::max(vtxBegin[1], vtxEnd[1]),
        std::max(vtxBegin[2], vtxEnd[2])
    };

    bcdef::PointRange pr{};

    /* Helper: set a varying axis' cell-center bounds. */
    auto set_var = [&](int axis, long long lo, long long hi, long long dim) {
        if (dim <= 1) {
            pr.begin[axis] = 1;
            pr.end[axis] = 1;
        } else {
            pr.begin[axis] = lo;
            pr.end[axis] = hi - 1;
        }
    };

    switch (face) {
        case bcdef::FaceDir::IMIN:
            pr.begin[0] = 1;
            pr.end[0] = 1;
            set_var(1, vmin[1], vmax[1], zone.nj());
            set_var(2, vmin[2], vmax[2], zone.nk());
            break;
        case bcdef::FaceDir::IMAX:
            pr.begin[0] = zone.ni() - 1;
            pr.end[0] = zone.ni() - 1;
            set_var(1, vmin[1], vmax[1], zone.nj());
            set_var(2, vmin[2], vmax[2], zone.nk());
            break;
        case bcdef::FaceDir::JMIN:
            pr.begin[1] = 1;
            pr.end[1] = 1;
            set_var(0, vmin[0], vmax[0], zone.ni());
            set_var(2, vmin[2], vmax[2], zone.nk());
            break;
        case bcdef::FaceDir::JMAX:
            pr.begin[1] = zone.nj() - 1;
            pr.end[1] = zone.nj() - 1;
            set_var(0, vmin[0], vmax[0], zone.ni());
            set_var(2, vmin[2], vmax[2], zone.nk());
            break;
        case bcdef::FaceDir::KMIN:
            pr.begin[2] = 1;
            pr.end[2] = 1;
            set_var(0, vmin[0], vmax[0], zone.ni());
            set_var(1, vmin[1], vmax[1], zone.nj());
            break;
        case bcdef::FaceDir::KMAX:
            pr.begin[2] = zone.nk() - 1;
            pr.end[2] = zone.nk() - 1;
            set_var(0, vmin[0], vmax[0], zone.ni());
            set_var(1, vmin[1], vmax[1], zone.nj());
            break;
    }

    return pr;
}

/*=====================================================================
  face_center_range (Plot3DZone overload)

  Same as the Zone overload, but uses Plot3DZone accessors.

  This is used by the Plot3D BC writer path to convert connectivity
  patches (in vertex space) into cell-center index ranges.
=====================================================================*/
bcdef::PointRange face_center_range(const bcdef::Plot3DZone& zone, bcdef::FaceDir face,
                                 const bcdef::IJK& vtxBegin, const bcdef::IJK& vtxEnd)
{
    bcdef::IJK vmin = {
        std::min(vtxBegin[0], vtxEnd[0]),
        std::min(vtxBegin[1], vtxEnd[1]),
        std::min(vtxBegin[2], vtxEnd[2])
    };
    bcdef::IJK vmax = {
        std::max(vtxBegin[0], vtxEnd[0]),
        std::max(vtxBegin[1], vtxEnd[1]),
        std::max(vtxBegin[2], vtxEnd[2])
    };

    bcdef::PointRange pr{};
    auto set_var = [&](int axis, long long lo, long long hi, long long dim) {
        if (dim <= 1) {
            pr.begin[axis] = 1;
            pr.end[axis] = 1;
        } else {
            pr.begin[axis] = lo;
            pr.end[axis] = hi - 1;
        }
    };

    switch (face) {
        case bcdef::FaceDir::IMIN:
            pr.begin[0] = 1;
            pr.end[0] = 1;
            set_var(1, vmin[1], vmax[1], zone.nj());
            set_var(2, vmin[2], vmax[2], zone.nk());
            break;
        case bcdef::FaceDir::IMAX:
            pr.begin[0] = zone.ni() - 1;
            pr.end[0] = zone.ni() - 1;
            set_var(1, vmin[1], vmax[1], zone.nj());
            set_var(2, vmin[2], vmax[2], zone.nk());
            break;
        case bcdef::FaceDir::JMIN:
            pr.begin[1] = 1;
            pr.end[1] = 1;
            set_var(0, vmin[0], vmax[0], zone.ni());
            set_var(2, vmin[2], vmax[2], zone.nk());
            break;
        case bcdef::FaceDir::JMAX:
            pr.begin[1] = zone.nj() - 1;
            pr.end[1] = zone.nj() - 1;
            set_var(0, vmin[0], vmax[0], zone.ni());
            set_var(2, vmin[2], vmax[2], zone.nk());
            break;
        case bcdef::FaceDir::KMIN:
            pr.begin[2] = 1;
            pr.end[2] = 1;
            set_var(0, vmin[0], vmax[0], zone.ni());
            set_var(1, vmin[1], vmax[1], zone.nj());
            break;
        case bcdef::FaceDir::KMAX:
            pr.begin[2] = zone.nk() - 1;
            pr.end[2] = zone.nk() - 1;
            set_var(0, vmin[0], vmax[0], zone.ni());
            set_var(1, vmin[1], vmax[1], zone.nj());
            break;
    }

    return pr;
}

} // namespace bcdef::boundary
