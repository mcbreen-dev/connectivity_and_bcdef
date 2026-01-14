# Testing

This project uses small shell scripts that drive `bc_define` and `bc_dump`
against a set of reference meshes.

## Full suite

```
./test/runtests_cgns.sh
```

- Writes output to `test/runtests.log` by default.
- Use `-v` to also echo output to the terminal.
- Each test copies a fresh mesh into `test/`, runs `bc_define`, then `bc_dump`.
- Temporary meshes are removed at the end of the run.

## Plot3D suite

```
./test/runtests_p3d.sh
```

- Writes output to `test/runtests.log` by default.
- Use `-v` to also echo output to the terminal.
- Each test copies a fresh mesh into `test/` and runs `bc_define`.
- Temporary meshes are removed at the end of the run.

## Single test

```
./test/run1test.sh <n>
```

Example:

```
./test/run1test.sh 3
```

- Writes output to `test/test<n>.log`.
- Copies a fresh mesh into `test/` and keeps it after the run.

## Suggested workflow

- Use `run1test.sh` while iterating on a single case.
- Run the full suite before tagging or sharing changes.
