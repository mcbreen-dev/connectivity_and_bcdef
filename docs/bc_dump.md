# bc_dump output

`bc_dump` prints a compact summary of zones, 1-to-1 connectivity, and BCs in a
CGNS file.

## Zones

The header lists the zone dimensions:

```
Zones in file "mesh.cgns":
Name      Ni   Nj   Nk
-------- ---- ---- ----
Zone1    51   51   1
```

## 1-to-1 connectivity

Connectivity entries show donor and receiver faces, receiver range, and the
transform vector:

```
1to1 links (2):
  • Donor Z3 i- → Receiver Z1 i+   r:[51,1,1]-[51,51,1]  T{1,2,0}
```

Notes:

- `r:[...]` is the receiver `PointRange`.
- For 2D zones, ranges are printed as `[i1,j1,1]-[i2,j2,1]`.
- `T{...}` is the CGNS transform vector.

## Boundary conditions

BC entries are printed per zone:

```
BCs (1):
  • inflow  [1,1,1]-[1,50,1]
```

Notes:

- The label uses the BC family name when present.
- For 2D zones, ranges are printed as `[i1,j1,1]-[i2,j2,1]`.
- When the BC is not `FamilySpecified`, the type is printed in parentheses.
