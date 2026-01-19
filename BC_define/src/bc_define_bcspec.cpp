/*─────────────────────────────────────────────────────────────
// File: src/bc_define_bcspec.cpp
// Boundary-condition specification parsing and validation
//─────────────────────────────────────────────────────────────
//
// This file implements parsing of user-provided boundary-condition
// definition files ("bcdef" files) and converts them into a structured,
// per-(zone,face) representation used by the BC_define pipeline.
//
// Responsibilities:
//   - Parse free-form, line-oriented BC specification text files.
//   - Normalize and resolve BC type tokens into CGNS BCType_t values.
//   - Interpret zone selectors and face selectors.
//   - Enforce consistency: a given (zone,face) may only be assigned
//     one BC type and family.
//   - Preserve stable CGNS FamilyName strings for later BC emission.
//
// Input format (conceptual):
//   Each non-empty, non-comment line describes a BC assignment:
//
//       <bc_type> <zone_selector> <face_selector>
//
//   Examples:
//       wall        1-4        i+
//       inflow      5,6        j-
//       symmetry    all        k+
//
// Design notes:
//   - Parsing is intentionally permissive: BC type tokens are matched
//     case-insensitively and may omit the "BC" prefix used by CGNS.
//   - Unknown BC tokens are mapped to BCTypeUserDefined rather than
//     rejected outright, allowing downstream policy decisions.
//   - Zone selectors are expanded eagerly so that conflicts can be
//     detected early and reported with precise diagnostics.
//
// Output:
//   parse_bcdef() returns an unordered_map keyed by (zone,face),
//   where each value records the CGNS BC type and the associated
//   FamilyName to be written later into the CGNS file.
//
// This module is purely concerned with *specification parsing*;
// it does not write CGNS nodes or interact with mesh connectivity.
───────────────────────────────────────────────────────────────*/
#include "bc_define_internal.hpp"

#include <fstream>
#include <sstream>

namespace fs::bc {

/*-------------------------------------------------------------
  parse_bc_type

  Convert a user-facing BC token into a CGNS BCType_t.

  This function is intentionally permissive:
    - It accepts common shorthand tokens (e.g. "wall", "symmetry")
    - It accepts common CGNS-style names with/without the "BC" prefix
      (e.g. "BCInflow" or "inflow")
    - It falls back to a brute-force scan over CGNS's known type names
      via cg_BCTypeName()

  Normalization:
    normalize_token() (from bc_define_internal.hpp) is used to make
    comparisons case-insensitive and whitespace/punctuation tolerant.

  Return value:
    - A concrete CGNS BCType_t for recognized tokens
    - BCTypeUserDefined if the token can't be resolved

  Note:
    Returning UserDefined (instead of throwing) keeps the parser robust
    and lets the caller decide whether unknown types are acceptable.
-------------------------------------------------------------*/
CGNS_ENUMT(BCType_t) parse_bc_type(const std::string& token)
{
    const std::string norm = normalize_token(token);

    /*---------------------------------------------------------
      Fast-path synonyms / shorthands

      These cover the most common BC labels used in typical bcdef files
      and allow convenient “plain English” inputs.
    ---------------------------------------------------------*/
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

    /*---------------------------------------------------------
      Slow-path: scan the CGNS BC type table.

      CGNS exposes a human-readable name for each enum value via
      cg_BCTypeName(). We normalize those names the same way and
      compare them to the user's token.

      We also accept tokens that omit the "BC" prefix:
        e.g. "WallInviscid" should match "BCWallInviscid".
    ---------------------------------------------------------*/
    for (int t = 0; t < NofValidBCTypes; ++t) {
        const char* name = cg_BCTypeName(static_cast<CGNS_ENUMT(BCType_t)>(t));
        if (!name)
            continue;
        std::string candidate = normalize_token(name);
        if (candidate == norm)
            return static_cast<CGNS_ENUMT(BCType_t)>(t);

        // If CGNS name begins with "bc...", allow the user to omit "bc".
        if (candidate.rfind("bc", 0) == 0 && candidate.substr(2) == norm)
            return static_cast<CGNS_ENUMT(BCType_t)>(t);
    }

    /*---------------------------------------------------------
      Unknown token: keep going, but mark as UserDefined.

      The higher-level pipeline can decide whether this should be:
        - accepted as BCTypeUserDefined
        - rejected with an error
        - mapped later via Family_t handling
    ---------------------------------------------------------*/
    return CGNS_ENUMV(BCTypeUserDefined);
}

/*-------------------------------------------------------------
  parse_bcdef

  Parse a BC definition file (bcdef) and produce a per-(zone,face)
  lookup table of desired BC assignments.

  File format (as implemented here):
    - Lines may include comments starting with '#'
    - After stripping comments and trimming:
        tokens[0]            : BC "type" token (e.g. wall, inflow)
        tokens[1..N-2] joined: zone selector/list token (no separators required)
        tokens[N-1]          : face selector token (e.g. i-, j+, k-)

  Example lines (conceptually):
      wall        1,2,3     i-
      inflow      4-8       j-
      farfield    all       i+

  The exact accepted zone selector grammar is implemented by
  parse_zone_list() in bc_define_internal.hpp.

  Validation / conflict handling:
    - If the same (zone,face) is specified multiple times, it must
      resolve to the same {family,type}. Otherwise an error is thrown.

  Family naming:
    - family_name is derived from the *first spelling* seen for a
      normalized BC token.
      Example: if the first file entry says "Wall", later "wall"
      entries will reuse "Wall" as the family name.
-------------------------------------------------------------*/
std::unordered_map<ZoneFaceKey, BCSpec, ZoneFaceKeyHash>
parse_bcdef(const std::string& path)
{
    std::ifstream in(path);
    if (!in)
        throw std::runtime_error("Could not open BC definition file: " + path);

    /*---------------------------------------------------------
      specs:
        Final output mapping (zone,face) -> BCSpec{family, type}

      family_names:
        Maps normalized bc_token -> canonical family name string.

      This ensures the same logical BC ("wall") produces a consistent
      FamilyName spelling even if the input varies in case.
    ---------------------------------------------------------*/
    std::unordered_map<ZoneFaceKey, BCSpec, ZoneFaceKeyHash> specs;
    std::unordered_map<std::string, std::string> family_names;

    std::string line;
    while (std::getline(in, line)) {

        /*-----------------------------------------------------
          Strip comments and ignore blank lines
        -----------------------------------------------------*/
        size_t hash = line.find('#');
        if (hash != std::string::npos)
            line = line.substr(0, hash);
        line = trim(line);
        if (line.empty())
            continue;

        /*-----------------------------------------------------
          Tokenize by whitespace

          Expected minimum:
            [bc_token] [zones_token...] [face_token]
          Anything less than 3 tokens is ignored (no-op).
        -----------------------------------------------------*/
        std::istringstream iss(line);
        std::vector<std::string> tokens;
        std::string tok;
        while (iss >> tok)
            tokens.push_back(tok);
        if (tokens.size() < 3)
            continue;

        /*-----------------------------------------------------
          Extract fields:
            - bc_token    : first token
            - face_token  : last token
            - zones_token : concatenation of middle tokens

          The concatenation allows users to write zone selectors with
          spaces (e.g., "1, 2, 3") without caring about whitespace.
        -----------------------------------------------------*/
        std::string bc_token = tokens.front();
        std::string face_token = tokens.back();
        std::string zones_token;
        for (size_t i = 1; i + 1 < tokens.size(); ++i)
            zones_token += tokens[i];

        /*-----------------------------------------------------
          Parse face selector and BC type
        -----------------------------------------------------*/
        fs::FaceDir face = fs::selector_to_facedir(face_token);
        CGNS_ENUMT(BCType_t) bc_type = parse_bc_type(bc_token);

        /*-----------------------------------------------------
          Canonicalize the family name for this BC token

          We key on the normalized bc token (case-insensitive, etc.)
          but store the original spelling from the first occurrence.
        -----------------------------------------------------*/
        std::string family_key = normalize_token(bc_token);
        if (family_names.find(family_key) == family_names.end())
            family_names[family_key] = bc_token;

        const std::string& family_name = family_names[family_key];

        /*-----------------------------------------------------
          Expand the zone selector into concrete zone indices,
          then populate the (zone,face)->BCSpec mapping.
        -----------------------------------------------------*/
        std::vector<int> zones = parse_zone_list(zones_token);
        for (int zone : zones) {
            ZoneFaceKey key{zone, face};

            // If this (zone,face) was already specified, it must agree.
            auto it = specs.find(key);
            if (it != specs.end()) {
                if (it->second.family != family_name || it->second.type != bc_type) {
                    throw std::runtime_error(
                        "Conflicting BCs for zone " + std::to_string(zone) + " face " +
                        fs::facedir_to_string(face));
                }
                continue;
            }

            // First time this (zone,face) is defined → store it.
            specs[key] = BCSpec{family_name, bc_type};
        }
    }

    return specs;
}

} // namespace fs::bc
