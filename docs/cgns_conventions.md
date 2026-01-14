# CGNS conventions used here

This project assumes a structured CGNS mesh with a single `Base_t` and one or
more `Zone_t` nodes.

## Dimensions

- `CellDimension` is read from `Base_t`.
- For 2D meshes, zones are treated as having `Nk = 1`.
- For 3D meshes, `Ni`, `Nj`, `Nk` are used as-is.

## Point ranges

- `PointRange` is the only point set type written for BCs and 1-to-1 links.
- For 2D meshes (index dimension = 2), CGNS stores 4 values for `PointRange`:
  `[i1, j1, i2, j2]`. Internally we track 6 values as
  `[i1, j1, 1]-[i2, j2, 1]` for consistency.
- For 3D meshes, 6 values are used: `[i1, j1, k1]-[i2, j2, k2]`.

## BC location

- BCs are written as `FamilySpecified` with `PointRange` and an explicit
  `GridLocation`.
- For 3D grids, `IFaceCenter/JFaceCenter/KFaceCenter` is used depending on face.
- For 2D grids, `EdgeCenter` is used.

## Families

- Each BC is associated with a `Family_t` node whose name is derived from the
  `bcdef.input` token (e.g., `inflow`).
- The BC nodes themselves are named `<family>_<face>_<count>` to keep them
  unique per zone+face.

## 1-to-1 connectivity

- Connectivity is written as `GridConnectivity1to1_t` with explicit receiver
  and donor `PointRange` values.
- A three-component transform vector is written to encode orientation.
- For 2D, the connectivity ranges are remapped to the 4-value CGNS form.
