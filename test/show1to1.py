#!/usr/bin/env python3
import os
import subprocess
import sys
from shutil import which


def _find_bc_dump(script_dir: str) -> str:
    candidates = [
        os.path.join(script_dir, "BC_define", "bin", "bc_dump"),
        os.path.join(script_dir, "..", "BC_define", "bin", "bc_dump"),
        os.path.join(script_dir, "..", "..", "BC_define", "bin", "bc_dump"),
    ]
    for path in candidates:
        if os.path.isfile(path):
            return os.path.abspath(path)
    path = which("bc_dump")
    if path:
        return path
    return os.path.abspath(candidates[0])


def main() -> int:
    if len(sys.argv) != 2:
        print("Usage: show1to1.py mesh.cgns")
        return 1

    mesh_path = sys.argv[1]
    if not os.path.isfile(mesh_path):
        print(f"Mesh not found: {mesh_path}")
        return 1

    script_dir = os.path.dirname(os.path.abspath(__file__))
    bc_dump = _find_bc_dump(script_dir)
    if not os.path.isfile(bc_dump):
        print(f"bc_dump not found at {bc_dump}. Build BC_define first.")
        return 1

    proc = subprocess.run([bc_dump, mesh_path], capture_output=True, text=True)
    output = proc.stdout.splitlines()
    if proc.returncode != 0:
        sys.stdout.write(proc.stdout)
        sys.stderr.write(proc.stderr)
        return proc.returncode

    keep = []
    for line in output:
        if line.startswith("Zones in file "):
            keep.append(line)
            continue
        if line.startswith("Name      Ni"):
            keep.append(line)
            continue
        if line.startswith("--------"):
            keep.append(line)
            continue
        if line.startswith("Zone"):
            keep.append(line)
            continue
        if "1to1 connections" in line:
            keep.append(line)
            continue
        if ("Donor " in line) and ("Receiver " in line):
            keep.append(line)
            continue

    print("\n".join(keep))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
