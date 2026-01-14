# bc_define algorithm

This document describes how `bc_define` detects connectivity and writes BCs.

## Overview

1. Open the mesh and read zone sizes.
2. Optionally purge existing `ZoneBC_t` and `ZoneGridConnectivity_t` nodes.
3. Build face-cell-center samples for each zone face.
4. Generate face-pair candidates with spatial hashing and bounds pruning.
5. Verify candidate pairs, solve orientation, and build 1-to-1 patches.
6. Record unmatched patches as boundary patches.
7. Apply BCs from `bcdef.input` (and optional autowall/autofarfield).
8. Optionally prune unused `Family_t` nodes (CGNS overwrite mode).

## Connectivity detection

The detector operates on face cell centers. For each structured zone:

### Face extraction
- For each exterior face, build a face-local grid of cell centers.
- For CGNS, the algorithm chooses between:
  - **Volume read**: read the full coordinate volume once, then compute face
    centers from the volume.
  - **Face-plane read**: read only the face planes.
- The choice uses `K2D`/`K3D` and `MEM_CAP_BYTES` from `include/k_values.hpp`.
  The heuristic compares volume element count against surface element count and
  also ensures the full read stays under a 1 GB cap.

### Candidate generation
- Each face has a bounding box and centroid computed from its cell centers.
- Faces are bucketed by their bounds (coarse spatial bins).
- Faces within the same bucket become candidate pairs.
- Candidates are pruned using face bounding-box overlap tests.

### Pair verification
For each candidate face pair:

1. Try to infer a UV mapping using a few matched samples.
2. If a valid mapping is found, use it to verify all points by direct index
   mapping (fast path).
3. If no mapping can be inferred, fall back to nearest-neighbor matching.

Pair verification is parallelized across face pairs. Outputs are sorted for
stable/deterministic results.

### Patch construction
- Matched points are grown into full-face or sub-face patches.
- Orientation is derived from index deltas to build a CGNS transform vector.
- Connectivity is recorded as explicit receiver/donor `PointRange` values.

### Left-handed checks
- For 3D blocks, a small local determinant is used to detect left-handed zones.
- Left-handed zones cause the run to fail.

## Boundary conditions

Boundary patches are assigned BCs in `bc_define`:

- The BC family and type are taken from `bcdef.input` for the matching
  zone+face. Families are created (or reused) under `Base_t` (CGNS).
- BCs are written as `FamilySpecified` + `PointRange`, with an explicit
  `GridLocation` for the face (CGNS).
- If `--autowall` is set, any boundary patch that still lacks a BC is assigned
  a `wall` family (`BCWall`).
- If `--autofarfield` is set, any boundary patch that still lacks a BC is
  assigned a `farfield` family (`BCFarfield`). These flags are mutually
  exclusive.
- Overlapping user-specified BCs on the same face are treated as an error.

## Non-overwrite mode (no duplicates)

When `--overwrite` is not set, the tool tries to avoid duplicates:

- Existing `GridConnectivity1to1_t` nodes are indexed by name and range. A new
  link is only written if no existing entry with the same name and range exists.
- Existing `BC_t` nodes are indexed by name and range. A new BC is only written
  if no existing entry with the same name and range exists.
- Autowall/autofarfield only applies to patches that do not already have any BC
  covering the same range.

## Overwrite mode cleanup

With `--overwrite`, the tool deletes `ZoneBC_t` and `ZoneGridConnectivity_t`
from each zone before writing. After writing, it removes any `Family_t` nodes
that are no longer referenced by zones or BCs.

## Plot3D outputs

For Plot3D meshes, no metadata can be embedded in the mesh, so `bc_define`
writes two files in the working directory:

- `1to1s`: 1-to-1 connectivity list.
- `bcs`: boundary condition list.

Format details are documented in `docs/plot3d_formats.md`.
