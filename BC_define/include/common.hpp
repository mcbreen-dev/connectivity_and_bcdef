/*─────────────────────────────────────────────────────────────
  File: include/common.hpp
  Lightweight types & helpers shared across BC_define modules
─────────────────────────────────────────────────────────────*/
#pragma once

#include <array>
#include <cstdint>
#include <string>
#include <unordered_map>
#include <stdexcept>

namespace fs {

/*----------- Face-direction enumeration --------------------*/
enum class FaceDir : uint8_t {
    IMIN, IMAX,
    JMIN, JMAX,
    KMIN, KMAX
};

/* Convert token like "i-" / "j+" / "k-" (case-insensitive) */
inline FaceDir selector_to_facedir(const std::string& tok)
{
    if (tok.size() != 2)
        throw std::invalid_argument("Selector must be exactly two chars (e.g. i-)");
    char axis = std::tolower(tok[0]);
    char pm   = tok[1];

    if (axis == 'i' && pm == '-') return FaceDir::IMIN;
    if (axis == 'i' && pm == '+') return FaceDir::IMAX;
    if (axis == 'j' && pm == '-') return FaceDir::JMIN;
    if (axis == 'j' && pm == '+') return FaceDir::JMAX;
    if (axis == 'k' && pm == '-') return FaceDir::KMIN;
    if (axis == 'k' && pm == '+') return FaceDir::KMAX;

    throw std::invalid_argument("Invalid selector token \"" + tok + '"');
}

/* Reverse mapping for pretty printing */
inline std::string facedir_to_string(FaceDir f)
{
    switch (f)
    {
        case FaceDir::IMIN: return "i-";
        case FaceDir::IMAX: return "i+";
        case FaceDir::JMIN: return "j-";
        case FaceDir::JMAX: return "j+";
        case FaceDir::KMIN: return "k-";
        case FaceDir::KMAX: return "k+";
        default:            return "??";
    }
}

/*----------- Index / range helpers --------------------------*/
using IJK = std::array<long long, 3>;   ///< 64-bit zone indices

struct PointRange           ///< CGNS-style PointRange(begin,end)
{
    IJK begin{};
    IJK end{};
};

/*----------- Geometry comparison tolerance ------------------*/
constexpr double GEOM_TOL = 1.0e-8;     ///< relative coordinate tolerence

/*----------- Error messages from CGNS -----------------------*/
#define CG_CALL(expr, msg)                                  \
    do {                                                    \
        if (int _cgstat = (expr)) {                         \
            cg_error_print();                               \
            throw std::runtime_error(msg);                  \
        }                                                   \
    } while (0)

} // namespace fs