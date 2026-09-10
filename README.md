# VASP Utility Scripts

Small utilities for VASP, molecular-dynamics, phonon, and NAMD calculations.

## Scripts

### `nve_band_analysis.py`

Analyze band energies from a continuous VASP MD `OUTCAR`.

Common usage:

```bash
python nve_band_analysis.py -b 313 316 317 321
python nve_band_analysis.py -b 316 317 -g 316 317
python nve_band_analysis.py -b 316 317 -g 316 317 --gap-fft
```

Here `-g A B` means the signed difference `E_B - E_A`. Use

```bash
python nve_band_analysis.py -h
```

for all options.

### `plot_band_population_checked.py`

Check and plot surfhop populations using the band basis stored in `HAMIL.h5`.

Example:

```bash
python plot_band_population_checked.py \
    -H HAMIL.h5 \
    -i excitation/averaged_results.h5 \
    -b 316 \
    -o excitation/band_population
```

Check without writing output files:

```bash
python plot_band_population_checked.py \
    -H HAMIL.h5 \
    -i excitation/averaged_results.h5 \
    -b 316 \
    -c
```

Default output is PNG + JSON. Add `--pdf` and/or `--csv` if needed.

Dependencies:

```text
numpy
h5py
matplotlib
```

### `xdatcar.py`

Split a VASP `XDATCAR` into individual POSCAR files:

```text
configs/POSCAR.0001
configs/POSCAR.0002
...
configs/POSCAR.3000
```

Set the number of frames with `NFRAME`.

Requirement:

```bash
pip install ase
```

Run with:

```bash
python3 xdatcar.py
```

### `init.sh`

Create one VASP run directory for each extracted structure:

```text
run/
├── 0001/
├── 0002/
├── ...
└── 3000/
```

Required files:

```text
INCAR
KPOINTS
POTCAR
configs/
```

Run with:

```bash
bash init.sh
```

`POSCAR` is copied into each directory. `KPOINTS` and `POTCAR` are symbolic links. `run/0001/INCAR` is a separate file with `ICHARG = 2`, while later frames link to the main `INCAR`, expected to use `ICHARG = 1`.

Typical workflow:

```bash
python3 xdatcar.py
bash init.sh
```

## Notes

- Band numbers use VASP's one-based band indices.
- Check input paths and calculation provenance before comparing different runs.
- Re-running `init.sh` does not remove existing VASP output files under `run/`.