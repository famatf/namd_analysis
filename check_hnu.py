import h5py
import numpy as np
from pathlib import Path

hnu = 1.5 # hnu in efield.rhai

with h5py.File("HAMIL.h5", "r") as f:
    bands = np.asarray(f["basis_list"][...]).reshape(-1)
    e = np.asarray(f["eig_t"][...])

    print("File:", Path("HAMIL.h5").resolve())
    print("Basis list:", bands.tolist())
    print("Energy array shape:", e.shape)
    print("scissor:", f["scissor"][...])

    assert e.ndim == 2 and e.shape[1] == len(bands), "Energy array shape mismatch with basis list"
    assert np.all(np.isfinite(e)), "Energy array contains non-finite values"

    print("\nBand energies (eV):")
    for i, band in enumerate(bands):
        print(f"Band {band}: mean={e[:,i].mean():.9f}",
              f"std={e[:,i].std():.9f}")

    print("\nEnergy differences (eV):")
    for a, b in [(316,317),(316,321),(317,321)]:
        if a not in bands or b not in bands:
            print(f"{a} -> {b}: not found in basis list")
            continue
        i = np.flatnonzero(bands == a)[0]
        j = np.flatnonzero(bands == b)[0]
        gap = e[:,j] - e[:,i]
        print(f"Band {a} -> {b}: mean={gap.mean():.9f}",
              f"std={gap.std():.9f}",
              f"min={gap.min():.9f}, max={gap.max():.9f}")
        if (a, b) == (316, 321):
            print(f"Photon energy: {hnu:.9f} eV")
            print(f"Photon energy - average energy: {hnu - gap.mean():+.9f} eV")