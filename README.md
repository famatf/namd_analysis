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

### `check_hnu.py`

Inspect band energies in `HAMIL.h5` and compare selected energy differences with a photon energy.

The script reports mean and standard deviation for every band in `basis_list`, then evaluates the fixed pairs `316 -> 317`, `316 -> 321`, and `317 -> 321` when they are present. For `316 -> 321`, it also prints the difference between the chosen photon energy and the average band-energy difference.

Set the photon energy directly in the script:

```
hnu = 1.5
```

Run with:

```
python check_hnu.py
```

The script expects `HAMIL.h5` in the current directory.

### `plot_gamma_phonon_projection.py`

Plot layer and polarization projections of Gamma-point phonon modes for a BP/MoS2 structure from Phonopy output.

The script reads eigenvectors from `qpoints.yaml` and uses atom symbols from that file or, if necessary, from a companion `phonopy.yaml`. It calculates BP/MoS2 layer weights and Cartesian polarization weights from normalized `|e|^2`.

Basic usage:

```
python plot_gamma_phonon_projection.py qpoints.yaml
```

For example, to label selected modes:

```
python plot_gamma_phonon_projection.py qpoints.yaml --label-modes 8,225,281
```

Default outputs are:

```
gamma_phonon_projection.png
gamma_phonon_projection.pdf
gamma_phonon_projection.csv
```

Dependencies:

```
numpy
matplotlib
PyYAML
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