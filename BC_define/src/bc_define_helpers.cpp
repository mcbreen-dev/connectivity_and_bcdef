/*─────────────────────────────────────────────────────────────
  File: src/bc_define_helpers.cpp
  Shared small utilities for bc_define
─────────────────────────────────────────────────────────────*/
#include "bc_define_internal.hpp"

#include <algorithm>
#include <cctype>
#include <filesystem>

namespace fs::bc {

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

std::string bc_patch_name(const std::string& family, fs::FaceDir face, int index)
{
    std::string name = family + "_" + fs::facedir_to_string(face) + "_" + std::to_string(index);
    if (name.size() > CGNS_MAX_NAME_LENGTH)
        name.resize(CGNS_MAX_NAME_LENGTH);
    return name;
}

fs::PointRange face_center_range(const fs::Zone& zone, fs::FaceDir face,
                                 const fs::IJK& vtxBegin, const fs::IJK& vtxEnd)
{
    fs::IJK vmin = {
        std::min(vtxBegin[0], vtxEnd[0]),
        std::min(vtxBegin[1], vtxEnd[1]),
        std::min(vtxBegin[2], vtxEnd[2])
    };
    fs::IJK vmax = {
        std::max(vtxBegin[0], vtxEnd[0]),
        std::max(vtxBegin[1], vtxEnd[1]),
        std::max(vtxBegin[2], vtxEnd[2])
    };

    fs::PointRange pr{};
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
        case fs::FaceDir::IMIN:
            pr.begin[0] = 1;
            pr.end[0] = 1;
            set_var(1, vmin[1], vmax[1], zone.nj());
            set_var(2, vmin[2], vmax[2], zone.nk());
            break;
        case fs::FaceDir::IMAX:
            pr.begin[0] = zone.ni() - 1;
            pr.end[0] = zone.ni() - 1;
            set_var(1, vmin[1], vmax[1], zone.nj());
            set_var(2, vmin[2], vmax[2], zone.nk());
            break;
        case fs::FaceDir::JMIN:
            pr.begin[1] = 1;
            pr.end[1] = 1;
            set_var(0, vmin[0], vmax[0], zone.ni());
            set_var(2, vmin[2], vmax[2], zone.nk());
            break;
        case fs::FaceDir::JMAX:
            pr.begin[1] = zone.nj() - 1;
            pr.end[1] = zone.nj() - 1;
            set_var(0, vmin[0], vmax[0], zone.ni());
            set_var(2, vmin[2], vmax[2], zone.nk());
            break;
        case fs::FaceDir::KMIN:
            pr.begin[2] = 1;
            pr.end[2] = 1;
            set_var(0, vmin[0], vmax[0], zone.ni());
            set_var(1, vmin[1], vmax[1], zone.nj());
            break;
        case fs::FaceDir::KMAX:
            pr.begin[2] = zone.nk() - 1;
            pr.end[2] = zone.nk() - 1;
            set_var(0, vmin[0], vmax[0], zone.ni());
            set_var(1, vmin[1], vmax[1], zone.nj());
            break;
    }

    return pr;
}

fs::PointRange face_center_range(const fs::Plot3DZone& zone, fs::FaceDir face,
                                 const fs::IJK& vtxBegin, const fs::IJK& vtxEnd)
{
    fs::IJK vmin = {
        std::min(vtxBegin[0], vtxEnd[0]),
        std::min(vtxBegin[1], vtxEnd[1]),
        std::min(vtxBegin[2], vtxEnd[2])
    };
    fs::IJK vmax = {
        std::max(vtxBegin[0], vtxEnd[0]),
        std::max(vtxBegin[1], vtxEnd[1]),
        std::max(vtxBegin[2], vtxEnd[2])
    };

    fs::PointRange pr{};
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
        case fs::FaceDir::IMIN:
            pr.begin[0] = 1;
            pr.end[0] = 1;
            set_var(1, vmin[1], vmax[1], zone.nj());
            set_var(2, vmin[2], vmax[2], zone.nk());
            break;
        case fs::FaceDir::IMAX:
            pr.begin[0] = zone.ni() - 1;
            pr.end[0] = zone.ni() - 1;
            set_var(1, vmin[1], vmax[1], zone.nj());
            set_var(2, vmin[2], vmax[2], zone.nk());
            break;
        case fs::FaceDir::JMIN:
            pr.begin[1] = 1;
            pr.end[1] = 1;
            set_var(0, vmin[0], vmax[0], zone.ni());
            set_var(2, vmin[2], vmax[2], zone.nk());
            break;
        case fs::FaceDir::JMAX:
            pr.begin[1] = zone.nj() - 1;
            pr.end[1] = zone.nj() - 1;
            set_var(0, vmin[0], vmax[0], zone.ni());
            set_var(2, vmin[2], vmax[2], zone.nk());
            break;
        case fs::FaceDir::KMIN:
            pr.begin[2] = 1;
            pr.end[2] = 1;
            set_var(0, vmin[0], vmax[0], zone.ni());
            set_var(1, vmin[1], vmax[1], zone.nj());
            break;
        case fs::FaceDir::KMAX:
            pr.begin[2] = zone.nk() - 1;
            pr.end[2] = zone.nk() - 1;
            set_var(0, vmin[0], vmax[0], zone.ni());
            set_var(1, vmin[1], vmax[1], zone.nj());
            break;
    }

    return pr;
}

} // namespace fs::bc
