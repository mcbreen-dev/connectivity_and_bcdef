/*─────────────────────────────────────────────────────────────
  File: include/common.hpp
  Types & helpers shared across BC_define modules

  This header defines:
  - Core enums and small structs used across the entire project
  - Lightweight geometry/index containers with fixed sizes
  - Face-direction encoding consistent with structured CGNS grids
  - A global geometric tolerance constant
  - A CGNS error-handling macro that converts C-style errors into
    C++ exceptions

  Design intent:
  - No heavy dependencies
  - No ownership of large data
  - Stable ABI-sized types suitable for hashing, copying, and logging

  Keeping these primitives here prevents cyclic dependencies and helps
  keep other headers focused on their own domains (mesh, connectivity, BCs).
─────────────────────────────────────────────────────────────*/
#pragma once

#include <array>
#include <cstdint>
#include <string>
#include <unordered_map>
#include <stdexcept>

namespace fs {

/*======================================================================
  FaceDir
======================================================================*/
/*
  Enumeration of the six boundary faces of a structured hexahedral zone (block).

    Each value represents a logical plane of constant index:
      IMIN / IMAX  : constant i index (left/right planes)
      JMIN / JMAX  : constant j index (front/back planes)
      KMIN / KMAX  : constant k index (bottom/top planes)

  Data type:
    - uint8_t (unsigned (non-negative) integer, 1 byte: 0-255)
        - Values are small, fixed, and frequently copied
        - uint8_t minimizes memory footprint when stored in large vectors
        (e.g., per-face metadata arrays)
        - Ordering is stable and relied upon in several switch statements

  This is used throughout:
  - Connectivity detection: describing faces being matched
  - BC definition: assigning BCs per zone face
  - Writers: choosing face-local parameter axes and transforms
*/
enum class FaceDir : uint8_t {
    IMIN, IMAX,
    JMIN, JMAX,
    KMIN, KMAX
};

/*======================================================================
  selector_to_facedir
======================================================================*/
/* Convert token like "i-" / "j+" / "k-" (case-insensitive) into a FaceDir*/
/*
  selector_to_facedir parses an *exactly two-character* selector:
    - First char: axis (i/j/k), case-insensitive
    - Second char: sign '-' for MIN, '+' for MAX

  Examples:
    "i-" => FaceDir::IMIN
    "J+" => FaceDir::JMAX

  Errors:
    - Throws std::invalid_argument if:
        * Token length ≠ 2
        * Axis letter is not i/j/k
        * Sign is not '+' or '-'

  This function is used when parsing bc-definition input
  files, where faces are specified in compact textual form.
*/
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

/*======================================================================
  facedir_to_string
======================================================================*/
/*
  Converts a FaceDir value back into its canonical two-character string
  representation ("i-", "j+", etc.).

  Return value:
    - std::string of length 2 for valid FaceDir values
    - "??" as a defensive fallback (should be unreachable)

  This function is used for:
    - Log messages
    - CGNS connectivity naming
    - Boundary-condition patch naming
*/
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

/*======================================================================
  Index / range helpers
======================================================================*/

/*----------- IJK ------------------*/
/*
  Fixed-size container for structured-grid indices (i, j, k).

  Type:
    std::array<long long, 3>

    - 3 × 8 bytes = 24 bytes

    - long long is used to safely store CGNS indices (cgsize_t-compatible)
    - Avoids overflow for large grids
    - Fixed-size array avoids heap allocation and allows trivial copying

  Indexing convention:
    - IJK[0] → i index
    - IJK[1] → j index
    - IJK[2] → k index

  Usage:
    - PointRange begin/end indices
    - Donor/receiver start indices in connectivity transforms
*/
using IJK = std::array<long long, 3>;

/*----------- PointRange ------------------*/
/*
  Represents a CGNS-style PointRange.

  Fields:
    begin : IJK
      - Lower corner of the range (inclusive), (i1, j1, k1)
    end   : IJK
      - Upper corner of the range (inclusive), (i2, j2, k2)

  Size:
    - 2 × IJK = 48 bytes total

  Notes:
    - Indices are expected to be 1-based when written to CGNS
    - In 2-D cases, the k-dimension is set to 1

  Usage:
    - Boundary-condition extents
    - GridConnectivity1to1 receiver/donor ranges
*/
struct PointRange
{
    IJK begin{};
    IJK end{};
};

/*----------- GEOM_TOL ------------------*/
/*
  Global geometric tolerance used in connectivity detection.

  Type:
    double (8 bytes)

  Default value:
    1.0e-8

  Usage:
    - Spatial hashing cell size selection
    - Squared-distance comparisons between face-center coordinates
    - Bounding-box overlap tests

    It is passed as fs::GEOM_TOL into ConnectivityDetector.

  Interpretation:
    - Absolute tolerance in coordinate space
    - Effectively assumes mesh coordinates are O(1)–O(10^3) in magnitude
*/
constexpr double GEOM_TOL = 1.0e-8;

/*----------- CG_CALL -----------------------*/
/*
  Macro wrapper for CGNS C API calls.

  Behavior:
    - Evaluates a CGNS function call
    - If the return value is non-zero:
        * Calls cg_error_print() to emit CGNS diagnostics
        * Throws std::runtime_error with the provided message

  Design notes:
    - Converts CGNS error handling into C++ exception flow
    - Avoids silent failures
    - Keeps calling code concise and readable

  Example:
    CG_CALL(cg_function(...), "message");

  Requirements:
    - This macro expects cg_error_print() to be available (from cgnslib.h).
*/
#define CG_CALL(expr, msg)                                  \
    do {                                                    \
        if (int _cgstat = (expr)) {                         \
            cg_error_print();                               \
            throw std::runtime_error(msg);                  \
        }                                                   \
    } while (0)

} // namespace fs