/*
  File: include/connectivity_writer.hpp

  Declarations for emitting CGNS GridConnectivity1to1_t nodes based on a detected
  structured 1:1 connection.

  Usage:
    - connectivity_core.cpp:
        - ConnectivityDetector::run(...) calls write_1to1(mesh, cp) for each
          ConnPatch produced by resolve_matches_cgns(...)

    - connectivity_verify.cpp:
        - resolve_matches_cgns(...) and resolve_matches_plot3d(...)
          construct ConnPatch records.
        - For CGNS output, resolve_matches_cgns sets:
            * recvZone/donorZone
            * recvDir/donorDir
            * offU/offV and sizeU/sizeV describing a patch in face-local UV
            * transform[0..2] as CGNS transform vector components
            * recvRange/donorRange as explicit PointRange (vertex indexing)
            * useExplicitRanges = true
          These ConnPatch objects are then written by connectivity_core.cpp.

    - connectivity_writer.cpp:
        - Implements write_1to1 and uses ConnPatch fields to compute:
            * receiver and donor PointRange arrays (cgsize_t[6])
            * transform vector (int[3])
          and then calls cg_1to1_write(...).

  Data type / size notes:
    - Zone indices: int (typically 32-bit). CGNS APIs take int zone indices.
    - FaceDir: enum class FaceDir : uint8_t (1 byte) from common.hpp.
    - offU/offV/sizeU/sizeV: long long (typically 64-bit), used for face-local
      indexing; later converted to cgsize_t when assembling cg_1to1_write ranges.
    - ori: uint8_t (1 byte) bitfield used only when useExplicitRanges==false
      (legacy/orientation-based writer path).
    - transform[3]: int[3] (typically 3 * 4 bytes). CGNS transform vector uses int.
    - PointRange indices: fs::IJK is std::array<long long,3> (3 * 8 bytes typical).
      These ranges are converted to cgsize_t arrays in connectivity_writer.cpp.
*/
#pragma once

#include "common.hpp"
#include <cstdint>

namespace fs {

class Mesh;

/*
  ConnPatch

  Represents a single structured 1:1 interface between:
    - a receiver zone face (recvZone, recvDir)
    - a donor zone face    (donorZone, donorDir)

  Two ways to describe the connectivity are supported by the writer:

    1) Explicit PointRanges (current path used by connectivity detection)
       - useExplicitRanges == true
       - recvRange and donorRange contain CGNS-style PointRange(begin,end)
         in index space compatible with cg_1to1_write.
       - transform[0..2] is already populated with the CGNS transform vector.
       - offU/offV/sizeU/sizeV/ori are not consulted by the writer in this mode.

    2) Face-local patch + orientation bits (fallback/legacy path)
       - useExplicitRanges == false
       - offU/offV/sizeU/sizeV describe a patch window in receiver face-local
         (U,V) coordinates; ori encodes swap/flip operations.
       - connectivity_writer.cpp computes receiver and donor PointRanges and
         transform from recvDir/donorDir + ori + patch window.

  Fields are populated as follows in this project:
    - resolve_matches_cgns(...) (src/connectivity_verify.cpp) sets:
        recvZone, donorZone, recvDir, donorDir,
        offU, offV, sizeU, sizeV,
        transform[0..2],
        recvRange, donorRange,
        useExplicitRanges = true
      and pushes ConnPatch into the list returned by ConnectivityDetector::run_collect.

    - resolve_matches_plot3d(...) also builds ConnPatch records, but those are
      written to text outputs (1to1s) by bc_define_plot3d.cpp; write_1to1 is not
      used for Plot3D output.
*/
struct ConnPatch
{
    /*
      recvZone / donorZone
        - 1-based zone indices (match fs::Zone::idx and CGNS zone numbering)
        - type: int
    */
    int      recvZone = 0;
    int      donorZone = 0;

    /*
      recvDir / donorDir
        - which face plane is connected on each zone
        - type: FaceDir (uint8_t underlying)
        - mapping between numeric face ids and FaceDir is in connectivity_utils.cpp
          via face_dir_from_id(...)
    */
    FaceDir  recvDir = FaceDir::IMIN;
    FaceDir  donorDir = FaceDir::IMIN;

    /*
      offU / offV
        - receiver-face-local starting offset (0-based) for a patch window
        - used only when useExplicitRanges == false
        - type: long long (64-bit typical)
    */
    long long offU = 0;
    long long offV = 0;

    /*
      sizeU / sizeV
        - receiver-face-local patch window size in vertex counts, not cell counts
          (connectivity_verify.cpp sets:
             sizeU = (uEnd - uStart + 2)
             sizeV = (vEnd - vStart + 2)
           where uEnd/uStart are cell indices on the face grid, so +2 converts
           a cell span into vertex span.)
        - used by writer when useExplicitRanges == false
        - type: long long (64-bit typical)
    */
    long long sizeU = 0;
    long long sizeV = 0;

    /*
      ori
        - 3-bit orientation encoding used by connectivity_writer.cpp when
          useExplicitRanges == false.

        Interpreted in connectivity_writer.cpp as:
          swapUV = ori & 4
          flipU  = ori & 1
          flipV  = ori & 2

        Type: uint8_t (1 byte).
    */
    uint8_t  ori = 0;

    /*
      transform
        - CGNS transform vector for cg_1to1_write:
            transform[rAxis] = ±(donorAxis+1) or 0
          where axis numbers are 1=I, 2=J, 3=K.
        - In the explicit-range path, resolve_matches_* computes and assigns it.
        - In the orientation-bit path, connectivity_writer.cpp derives it.
        - Stored as int[3] (typically 12 bytes).
    */
    int      transform[3] = {0, 0, 0};

    /*
      recvRange / donorRange
        - CGNS-style PointRange(begin,end) as fs::PointRange
        - Each endpoint is fs::IJK = std::array<long long,3>
        - These are written by connectivity_writer.cpp when useExplicitRanges==true.
    */
    PointRange recvRange{};
    PointRange donorRange{};

    /*
      useExplicitRanges
        - If true:
            - writer uses recvRange/donorRange and transform directly
            - offU/offV/sizeU/sizeV/ori are ignored
        - If false:
            - writer computes ranges/transform from face-local patch description
    */
    bool    useExplicitRanges = false;
};

/*
  write_1to1(mesh, cp)

  Writes one CGNS GridConnectivity1to1_t under:
    Base_t / Zone_t(cp.recvZone) / ZoneGridConnectivity_t

  Implementation (src/connectivity_writer.cpp):
    - Computes receiver and donor PointRange arrays (cgsize_t[6]) and transform (int[3])
      either from explicit ranges or from face-local patch + ori.
    - Reads existing 1to1 connections once per receiver zone to avoid duplicates.
    - Generates a unique connection name based on:
        "Z<recv>_<recvDir>_to_Z<donor>_<donorDir>[_<n>]"
    - Calls cg_1to1_write(...). Throws std::runtime_error on CGNS errors.
*/
void write_1to1(Mesh& mesh, const ConnPatch& cp);

} // namespace fs
