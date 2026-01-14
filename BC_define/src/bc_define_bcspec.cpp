/*─────────────────────────────────────────────────────────────
  File: src/bc_define_bcspec.cpp
  BC definition parsing and validation
─────────────────────────────────────────────────────────────*/
#include "bc_define_internal.hpp"

#include <fstream>
#include <sstream>

namespace fs::bc {

CGNS_ENUMT(BCType_t) parse_bc_type(const std::string& token)
{
    const std::string norm = normalize_token(token);
    if (norm == "inflow" || norm == "bcinflow")
        return CGNS_ENUMV(BCInflow);
    if (norm == "outflow" || norm == "bcoutflow")
        return CGNS_ENUMV(BCOutflow);
    if (norm == "zerograd" || norm == "zerogradient")
        return CGNS_ENUMV(BCExtrapolate);
    if (norm == "symmetry" || norm == "symmetryplane")
        return CGNS_ENUMV(BCSymmetryPlane);
    if (norm == "wall" || norm == "bcwall")
        return CGNS_ENUMV(BCWall);
    if (norm == "wallviscous")
        return CGNS_ENUMV(BCWallViscous);
    if (norm == "wallinviscid")
        return CGNS_ENUMV(BCWallInviscid);
    if (norm == "farfield")
        return CGNS_ENUMV(BCFarfield);
    if (norm == "extrapolate")
        return CGNS_ENUMV(BCExtrapolate);
    if (norm == "neumann")
        return CGNS_ENUMV(BCNeumann);
    if (norm == "dirichlet")
        return CGNS_ENUMV(BCDirichlet);
    if (norm == "general")
        return CGNS_ENUMV(BCGeneral);
    if (norm == "userdefined")
        return CGNS_ENUMV(BCTypeUserDefined);

    for (int t = 0; t < NofValidBCTypes; ++t) {
        const char* name = cg_BCTypeName(static_cast<CGNS_ENUMT(BCType_t)>(t));
        if (!name)
            continue;
        std::string candidate = normalize_token(name);
        if (candidate == norm)
            return static_cast<CGNS_ENUMT(BCType_t)>(t);
        if (candidate.rfind("bc", 0) == 0 && candidate.substr(2) == norm)
            return static_cast<CGNS_ENUMT(BCType_t)>(t);
    }

    return CGNS_ENUMV(BCTypeUserDefined);
}

std::unordered_map<ZoneFaceKey, BCSpec, ZoneFaceKeyHash>
parse_bcdef(const std::string& path)
{
    std::ifstream in(path);
    if (!in)
        throw std::runtime_error("Could not open BC definition file: " + path);

    std::unordered_map<ZoneFaceKey, BCSpec, ZoneFaceKeyHash> specs;
    std::unordered_map<std::string, std::string> family_names;

    std::string line;
    while (std::getline(in, line)) {
        size_t hash = line.find('#');
        if (hash != std::string::npos)
            line = line.substr(0, hash);
        line = trim(line);
        if (line.empty())
            continue;

        std::istringstream iss(line);
        std::vector<std::string> tokens;
        std::string tok;
        while (iss >> tok)
            tokens.push_back(tok);
        if (tokens.size() < 3)
            continue;

        std::string bc_token = tokens.front();
        std::string face_token = tokens.back();
        std::string zones_token;
        for (size_t i = 1; i + 1 < tokens.size(); ++i)
            zones_token += tokens[i];

        fs::FaceDir face = fs::selector_to_facedir(face_token);
        CGNS_ENUMT(BCType_t) bc_type = parse_bc_type(bc_token);

        std::string family_key = normalize_token(bc_token);
        if (family_names.find(family_key) == family_names.end())
            family_names[family_key] = bc_token;

        const std::string& family_name = family_names[family_key];
        std::vector<int> zones = parse_zone_list(zones_token);
        for (int zone : zones) {
            ZoneFaceKey key{zone, face};
            auto it = specs.find(key);
            if (it != specs.end()) {
                if (it->second.family != family_name || it->second.type != bc_type) {
                    throw std::runtime_error(
                        "Conflicting BCs for zone " + std::to_string(zone) + " face " +
                        fs::facedir_to_string(face));
                }
                continue;
            }
            specs[key] = BCSpec{family_name, bc_type};
        }
    }

    return specs;
}

} // namespace fs::bc
