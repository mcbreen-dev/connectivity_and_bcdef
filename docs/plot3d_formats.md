# Plot3D output formats

When `bc_define` reads a Plot3D `.x` mesh, it writes two files in the current
working directory: `1to1s` and `bcs`.

These files are intended to be easy to parse while remaining human-readable.
All indices are 1-based and inclusive.

## 1to1s

Header:
- Comments starting with `#` include zone sizes and column definitions.

Columns:
```
recv_zone donor_zone
r_i1 r_j1 r_k1 r_i2 r_j2 r_k2
d_i1 d_j1 d_k1 d_i2 d_j2 d_k2
t1 t2 t3
```

Where:
- `recv_zone` is the receiver zone index.
- `donor_zone` is the donor zone index.
- `r_*` and `d_*` are receiver/donor point ranges.
- `t1 t2 t3` is the transform vector (same semantics as CGNS).

## bcs

Header:
- Comments starting with `#` include zone sizes and column definitions.

Columns:
```
zone family i1 j1 k1 i2 j2 k2
```

Where:
- `zone` is the zone index.
- `family` is the BC family name (string).
- `i1 j1 k1 i2 j2 k2` is the boundary point range.
