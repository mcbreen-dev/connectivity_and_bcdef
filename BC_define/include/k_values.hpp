/*
  File: include/k_values.hpp

  Heuristic constants for choosing between two CGNS coordinate read strategies
  in src/connectivity_face_build.cpp (build_faces_cgns):

    Strategy A: read 6 boundary face planes per zone (FacePlane)
      - reads CoordinateX/Y/Z for one constant-index plane at a time

    Strategy B: read full zone volume once (VolumeCoords)
      - reads CoordinateX/Y/Z for the full (ni*nj*nk) vertex volume
      - face centers are then computed by indexing into the volume arrays

  build_faces_cgns selects Strategy B when both conditions hold for the zone:

      vol_elems  <=  k * face_elems
      vol_bytes  <=  MEM_CAP_BYTES

  where:
    - vol_elems  = ni * nj * nk
    - face_elems = sum of boundary plane vertex counts used for center evaluation
        * 2D (nk==1): 2*(ni + nj)
        * 3D:         2*(ni*nj + ni*nk + nj*nk)
    - k is K2D when nk==1, otherwise K3D
    - vol_bytes = vol_elems * 3 * sizeof(double)
        (three coordinate arrays, double precision)

  These constants are used only by:
    - src/connectivity_face_build.cpp

  Data types:
    - K2D, K3D are double because they participate in floating-point comparisons.
    - MEM_CAP_BYTES is unsigned long long to hold 1 GiB exactly without overflow.

  Sizes:
    - MEM_CAP_BYTES = 1ULL << 30 = 1,073,741,824 bytes (1 GiB).
*/
#pragma once

namespace fs {

/*
  K2D

  Threshold multiplier used when zone nk==1 (2D structured zones).

  Applied in:
    use_volume = (vol_elems <= K2D * face_elems) && (vol_bytes <= MEM_CAP_BYTES)
*/
constexpr double K2D = 2.681; //Default found by scripts/compute_k_values.sh on my 2025 MBP (128GB memory, M4 Max w/ 16 CPU cores)

/*
  K3D

  Threshold multiplier used when zone nk>1 (3D structured zones).

  Applied in:
    use_volume = (vol_elems <= K3D * face_elems) && (vol_bytes <= MEM_CAP_BYTES)
*/
constexpr double K3D = 60.253; //Default found by scripts/compute_k_values.sh on my 2025 MBP (128GB memory, M4 Max w/ 16 CPU cores)

/*
  MEM_CAP_BYTES

  Maximum allowed memory footprint for VolumeCoords in build_faces_cgns.

  This is compared against:
    vol_bytes = (ni*nj*nk) * 3 * sizeof(double)

  With sizeof(double)=8, the maximum volume element count allowed by this cap is:
    MEM_CAP_BYTES / (3*8) = 1,073,741,824 / 24 = 44,739,242 elements (floor).
*/
constexpr unsigned long long MEM_CAP_BYTES = 1ULL << 30;
}
