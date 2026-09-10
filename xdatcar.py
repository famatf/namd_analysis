#!/usr/bin/env python3
"""
Script to convert XDATCAR to a series of POSCAR files where each POSCAR corresponds to a frame in the XDATCAR file.
The number of frames to extract can be set by changing the NFRAME variable.
The output POSCAR files will be saved in a directory named 'configs'.
"""
from pathlib import Path
from ase.io import read, write

XDATCAR = "XDATCAR"
NFRAME = 3000               # Should be set as you want, but not more than the number of frames in XDATCAR
OUTDIR = Path("configs")

print(f"Reading {XDATCAR} ...")
configs = read(XDATCAR, format="vasp-xdatcar", index=":")

print(f"Frames found: {len(configs)}")

if len(configs) < NFRAME:
    raise RuntimeError(
        f"XDATCAR only contains {len(configs)} frames, "
        f"but {NFRAME} frames are required."
    )

OUTDIR.mkdir(exist_ok=True)

for i, atoms in enumerate(configs[:NFRAME], start=1):
    outfile = OUTDIR / f"POSCAR.{i:04d}"
    write(
        outfile,
        atoms,
        format="vasp",
        vasp5=True,
        direct=True
    )

    if i % 100 == 0 or i == 1 or i == NFRAME:
        print(f"Wrote {outfile}")

print(f"Done: {NFRAME} configurations written to {OUTDIR}/")