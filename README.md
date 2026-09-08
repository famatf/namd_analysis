# namd_analysis scripts
Analysis and plotting scripts for calculations

## Scripts Usage (selected)

### plot_band_population_checked.py

#### Input

Use one continuous MD `OUTCAR` containing `E-fermi` records and band-energy tables for the requested k-point and spin component. The default input is `OUTCAR` in the current working directory.

```text
nve/
├── OUTCAR
└── nve_band_analysis.py
```

Band numbers refer to the one-based indices printed by VASP. The numbers used in the examples are illustrative; replace them with those relevant to your calculation. The script does not assign material, layer, or orbital character to a band.

#### Plot selected bands

```bash
python nve_band_analysis.py -b 313 316 317 321
```

This saves `nve_bands_energies.png`. FFT is disabled by default.

#### Calculate individual-band FFTs

```bash
python nve_band_analysis.py -b 313 316 317 321 --fft
```

This also saves `nve_bands_band_fft.png`, with the selected bands on a shared spectrum plot. `-f` is an alias for `--fft`.

#### Analyze an energy difference without FFT

```bash
python nve_band_analysis.py -b 316 317 -g 316 317
```

The order matters: **`-g A B` means `E_B - E_A`**, not the absolute difference. This example plots `E317 - E316` and saves `nve_bands_gap.png` in addition to the band-energy plot.

#### Calculate the energy-difference FFT

```bash
python nve_band_analysis.py -b 316 317 -g 316 317 --gap-fft
```

This produces the band-energy trace, the energy-difference trace, and its FFT. Individual-band FFT remains disabled unless `--fft` is also supplied.

The two indices supplied to `-g` are automatically included when reading the input, even if they are absent from `-b`. The band-energy plot and individual-band FFT still use only the bands selected with `-b`.

#### Enable both FFT types and save every format

```bash
python nve_band_analysis.py -i OUTCAR -b 313 316 317 321 -f -g 316 317 --gap-fft --save png pdf csv json -o results/nve
```

The output directory is created if needed.

#### Export numerical data only

```bash
python nve_band_analysis.py -b 316 317 -g 316 317 --gap-fft --save csv json -o results/nve
```

`--save` replaces the format selection. For example, `--save pdf` produces PDF plots without PNG, CSV, or JSON files. The default is **PNG only**.

#### Select a frame range and frequency display range

```bash
python nve_band_analysis.py -b 313 316 317 321 -f --start 501 --end 2500 --fmax 600
```

Frame bounds are one-based and inclusive. `--fmax` limits only the FFT plot display, in cm⁻¹; exported FFT data retain the full frequency range.

#### Command-line options

```bash
python nve_band_analysis.py -h
```

| Option                           | Description                                        | Default                     |
| -------------------------------- | -------------------------------------------------- | --------------------------- |
| `-i`, `--outcar PATH`            | Input OUTCAR path                                  | `OUTCAR`                    |
| `-b`, `--bands B1 B2 ...`        | Unique positive band indices, separated by spaces  | Required                    |
| `-f`, `--fft`                    | Calculate an FFT for each band selected with `-b`  | Off                         |
| `-g`, `--gap A B`                | Analyze the signed difference `E_B - E_A`          | Off                         |
| `--gap-fft`                      | Calculate the energy-difference FFT; requires `-g` | Off                         |
| `--save`, `--formats FORMAT ...` | Select `png`, `pdf`, `csv`, and/or `json`          | `png`                       |
| `-o`, `--prefix PREFIX`          | Output prefix, optionally including a directory    | `nve_bands`                 |
| `--start N`                      | First retained frame, one-based and inclusive      | `1`                         |
| `--end N`                        | Last retained frame, inclusive                     | Last available frame        |
| `--dt FS`                        | Uniform saved-frame spacing in fs; overrides POTIM | First POTIM found in OUTCAR |
| `-k`, `--kpoint N`               | Selected k-point index                             | `1`                         |
| `-s`, `--spin N`                 | Selected spin-component index                      | `1`                         |
| `--reference MODE`               | `mean-fermi` or `raw`                              | `mean-fermi`                |
| `--window MODE`                  | FFT window: `hann` or `none`                       | `hann`                      |
| `--detrend MODE`                 | FFT preprocessing: `mean` or `linear`              | `mean`                      |
| `--fmax VALUE`                   | Upper FFT plot limit in cm⁻¹                       | Full computed range         |
| `--dpi N`                        | Resolution passed to plot export                   | `400`                       |
| `-h`, `--help`                   | Display help and examples                          | —                           |

#### Output files

For `-o results/nve`, filenames follow this scheme. Each format is written only when selected with `--save`.

| Filename                             | Content                                                      | Required analysis option |
| ------------------------------------ | ------------------------------------------------------------ | ------------------------ |
| `nve_energies.png` / `.pdf`          | Selected band energies on one shared axis                    | Always                   |
| `nve_energies.csv`                   | Frame indices, times, Fermi energies, raw/referenced band energies, and occupations | Always                   |
| `nve_gap.png` / `.pdf` / `.csv`      | Signed energy difference versus time                         | `-g A B`                 |
| `nve_band_fft.png` / `.pdf` / `.csv` | Individual-band amplitude spectra                            | `--fft`                  |
| `nve_gap_fft.png` / `.pdf` / `.csv`  | Energy-difference amplitude spectrum                         | `--gap-fft`              |
| `nve_analysis.json`                  | Input path, settings, frame selection, reference, FFT metadata, and output paths | Always                   |

The energy CSV includes any additional bands needed for the requested energy difference. Missing occupation values are represented as `nan`. FFT CSV files contain frequency in THz, wavenumber in cm⁻¹, and amplitudes in eV.

JSON stores metadata and, when applicable, energy-difference statistics: mean, population standard deviation (`ddof=0`), minimum, maximum, and peak-to-peak range. Full time-series and spectral arrays are exported through CSV, not JSON.

Existing files with the same output names are overwritten. Files from earlier runs that are not selected in the current run are not removed.

#### Energy and time conventions

##### Energy reference

By default, one constant—the mean Fermi energy over the selected frames—is subtracted from every band:

```text
E_plot,n(t) = E_n(t) - mean(E_F)
```

This is not a frame-by-frame subtraction of `E_F(t)`, nor is each energy trace centered independently for plotting. Use `--reference raw` to plot the stored energies without this subtraction.

Energy differences are calculated directly from the raw energies, so the common reference cancels. No scissor correction, smoothing, or interpolation is applied.

##### Time axis

The first selected retained frame is assigned `0 fs`. For `N` selected frames separated by `dt`, the last plotted time is `(N - 1) × dt`.

The script uses `--dt` when provided; otherwise, it uses the first `POTIM` found in the file. `--dt` is the spacing between **saved energy frames**. If energies were saved every several MD steps, supply the corresponding saved-frame spacing. This does not repair missing or irregularly spaced frames.

#### FFT method

Individual-band and energy-difference FFTs are independent options. FFT requires at least four selected frames; a time-series plot requires at least two.

For each requested series, the script:

1. Removes its mean.
2. Optionally removes a least-squares linear trend with `--detrend linear`.
3. Applies a Hann window by default, or a rectangular window with `--window none`.
4. Computes a real FFT using `numpy.fft.rfft`.
5. Divides the magnitude by the window sum and doubles the positive-frequency bins, except the Nyquist bin when present.

The result is a **one-sided amplitude spectrum in eV**, not a power spectral density. Curves are not individually rescaled to a maximum of one. The default mean removal acts before windowing; a nonzero DC bin can still remain after applying the window. Amplitudes of off-bin components can depend on window choice and spectral leakage.

Frequency conversion uses the speed of light `c = 2.99792458 × 10^10 cm/s`:

```text
frequency_THz = frequency_cycles_per_fs × 1000
wavenumber_cm^-1 = frequency_cycles_per_fs × 10^15 / c
```

FFT bin spacing is `10^15 / (N × dt × c)` cm⁻¹. The theoretical Nyquist limit is `10^15 / (2 × dt × c)` cm⁻¹. No zero-padding is applied. Bin spacing is not a peak-position uncertainty estimate, and windowing affects spectral resolution.

An FFT peak describes a periodic component of the selected energy signal. The script does not identify phonon modes or assign mode symmetry from a peak.

#### Parsing assumptions and limitations

- The input is one continuous MD OUTCAR, not a collection of snapshot directories or an automatically joined set of restarted runs.
- A target band table is associated with the most recent `E-fermi` record. If several complete target tables share that record, the last complete table is retained.
- The time axis assumes one retained final eigenvalue table per MD frame. Table completeness alone does not establish SCF convergence, verify the NVE ensemble, or detect completely omitted frames.
- An incomplete trailing target table is ignored with a warning. An incomplete target table before a later complete table causes an error rather than silently compressing the time axis.
- The parser requires the requested band indices to be present. When `NBANDS` is available, the table must contain all indices from `1` through `NBANDS`.
- The selected k-point and spin component are fixed throughout the analysis. Compatibility with every VASP output variant, including all spinor/noncollinear layouts, has not been established.
- Band indices are not tracked by wavefunction overlap or orbital character. Band crossings can change the physical character associated with a fixed index.
- The file must contain `E-fermi` records even when `--reference raw` is selected, because the parser uses those records to group tables.

#### Troubleshooting

| Message or symptom                            | What to check                                                |
| --------------------------------------------- | ------------------------------------------------------------ |
| `Cannot find NVE OUTCAR`                      | Run in the input directory or specify `-i /path/to/OUTCAR`.  |
| `No complete target band tables were found`   | Confirm that eigenvalues are present and that band, k-point, spin, and NBANDS values match the file. |
| Incomplete target table inside the trajectory | Inspect the affected input section; the script refuses to skip an interior incomplete table. |
| Missing or invalid POTIM                      | Supply a verified, uniform saved-frame spacing with `--dt`.  |
| `--gap-fft requires -g A B`                   | Specify the pair, for example `-g 316 317 --gap-fft`.        |
| No PDF, CSV, or JSON was created              | These formats are opt-in; select them with `--save`.         |
| FFT plot extends to very high wavenumbers     | Use `--fmax`, for example `--fmax 600`, to restrict the displayed range. |

#### Validation

The current implementation was checked with synthetic OUTCAR data for energy extraction, energy-difference cancellation, FFT frequency conversion and amplitude scaling, Hann-window correction, even/odd sample counts, the Nyquist bin, incomplete tables, frame selection, independent FFT options, and output-format selection. These checks do not replace verification against the structure and sampling of your own OUTCAR.
