#!/usr/bin/env python3
"""Check/plot surfhop populations using an explicitly selected HAMIL basis.

From the calculation directory containing HAMIL.h5:
  python plot_band_population_checked.py --hamil HAMIL.h5 --input excitation/averaged_results.h5 --expected-initial-band 316 --check-only

By default, plotting writes PNG and provenance JSON.
Add --pdf and/or --csv to export those optional formats.
HDF5 files are opened read-only. Matching dimensions alone do not prove that
the selected HAMIL generated the selected result; check the run configuration.
Dependencies: numpy, h5py; matplotlib is needed only for plotting.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import h5py
import numpy as np


def find_dataset(handle, name, optional=False):
    if name in handle and isinstance(handle[name], h5py.Dataset):
        return handle[name]
    matches = []
    handle.visititems(lambda path, obj: matches.append(path)
                      if isinstance(obj, h5py.Dataset)
                      and path.rsplit('/', 1)[-1] == name else None)
    if len(matches) == 1:
        return handle[matches[0]]
    if not matches and optional:
        return None
    raise ValueError('Dataset {!r}: expected one match, found {}'.format(name, matches))


def read_basis(values):
    values = np.asarray(values)
    if (values.ndim != 1 or values.size == 0
            or values.dtype.kind not in 'iu'
            or np.any(values <= 0)
            or np.unique(values).size != values.size):
        raise ValueError('basis_list must be a nonempty vector of unique positive integer band IDs.')
    return values.astype(int)


def orient_population(values, nbands, time_axis=None):
    values = np.asarray(values)
    if values.ndim > 2:
        values = np.squeeze(values)
    if values.ndim != 2 or 0 in values.shape:
        raise ValueError('Expected a nonempty 2-D population array; got {}'.format(values.shape))
    if np.iscomplexobj(values):
        if not np.all(np.isfinite(values)) or np.max(np.abs(values.imag)) > 1e-10:
            raise ValueError('Populations contain invalid values or a significant imaginary component.')
        values = values.real
    if not np.all(np.isfinite(values)):
        raise ValueError('Populations contain NaN or infinity.')
    if time_axis is None:
        band_axes = [axis for axis, length in enumerate(values.shape) if length == nbands]
        if len(band_axes) != 1:
            raise ValueError('HAMIL has {} bands but population shape is {}. '
                             'No unique band axis: check result provenance; for a square array, '
                             'specify --time-axis from the data definition. '
                             'No columns have been discarded or relabeled.'.format(nbands, values.shape))
        time_axis = 1 - band_axes[0]
    if values.shape[1 - time_axis] != nbands:
        raise ValueError('The selected population band axis does not match HAMIL basis_list.')
    return np.asarray(values if time_axis == 0 else values.T, dtype=float), time_axis


def read_time(handle, count, dataset_name=None, dt=None):
    names = [dataset_name] if dataset_name else ['time_fs', 'time', 'times', 't']
    candidates = []
    for name in names:
        ds = find_dataset(handle, name, optional=dataset_name is None)
        if ds is not None:
            candidates.append(ds)
    if len(candidates) > 1:
        raise ValueError('Multiple time datasets found; select one using --time-dataset.')
    if candidates:
        ds = candidates[0]
        raw = np.asarray(ds[...]).squeeze()
        if np.iscomplexobj(raw) or raw.ndim != 1 or raw.size != count:
            raise ValueError('Time dataset {} shape {} does not match {} samples.'.format(ds.name, raw.shape, count))
        times = np.asarray(raw, dtype=float)
        unit = ds.attrs.get('units', ds.attrs.get('unit'))
        if isinstance(unit, bytes):
            unit = unit.decode('utf-8')
        if unit is not None and str(unit).strip().lower() not in ['fs', 'femtosecond', 'femtoseconds']:
            raise ValueError('Time dataset units {!r} are not fs; confirm the conversion before plotting.'.format(unit))
        source = '{} ({})'.format(ds.name, 'units=fs' if unit is not None else 'fs assumed; no unit attribute')
    elif dt is not None:
        times = np.arange(count, dtype=float) * dt
        source = 'generated from explicitly supplied dt={} fs; origin=0'.format(dt)
    else:
        return None, 'No time dataset found; sample indices only. Plotting requires a verified --dt in fs.'
    if not np.all(np.isfinite(times)) or np.any(np.diff(times) <= 0):
        raise ValueError('Time must be finite and strictly increasing.')
    return times, source


def decode_efield(values):
    values = np.asarray(values)
    if values.dtype == np.uint8:
        return values.tobytes().decode('utf-8', errors='replace')
    if values.ndim == 0:
        value = values.item()
        return value.decode('utf-8', errors='replace') if isinstance(value, bytes) else str(value)
    return str(values.tolist())


def plot_population(times, bands, pops, stem, dpi, save_pdf=False):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 11})
    fig, ax = plt.subplots(figsize=(8.2, 5.0))
    # Index discrete palettes directly: never stretch ten colors over sixteen bands.
    if len(bands) <= 10:
        colors = plt.get_cmap('tab10').colors
    elif len(bands) <= 20:
        colors = plt.get_cmap('tab20').colors
    else:
        colors = plt.get_cmap('hsv')(np.arange(len(bands)) / len(bands))
    for col, band in enumerate(bands):
        main = int(band) in (316, 317)
        ax.plot(times, pops[:, col], color=colors[col], label='Band {}'.format(band),
                linewidth=2.4 if main else 1.25, alpha=1 if main else .8,
                zorder=3 if main else 2)
    ax.set_xlabel('Time (fs)')
    ax.set_ylabel('Population')
    if times.size > 1:
        ax.set_xlim(times[0], times[-1])
    ax.set_ylim(min(0, pops.min()) - .03, max(1, pops.max()) + .03)
    ax.grid(True, linestyle='--', linewidth=.55, alpha=.32)
    ax.legend(loc='center left', bbox_to_anchor=(1.01, .5), frameon=False,
              ncol=1 if len(bands) <= 12 else 2, fontsize=9)
    fig.tight_layout()
    for suffix in (('.png', '.pdf') if save_pdf else ('.png',)):
        fig.savefig(str(stem) + suffix, dpi=dpi, bbox_inches='tight')
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""Examples (run from the directory containing HAMIL.h5):
  Short form:
    python plot_band_population_checked.py -H HAMIL.h5 -i excitation/averaged_results.h5 -b 316 -o excitation/band_population

  Full form (existing commands still work):
    python plot_band_population_checked.py --hamil HAMIL.h5 --input excitation/averaged_results.h5 --expected-initial-band 316 --output excitation/band_population

  Check only, without writing files:
    python plot_band_population_checked.py -H HAMIL.h5 -i excitation/averaged_results.h5 -b 316 -c

  Also save PDF and CSV:
    python plot_band_population_checked.py -H HAMIL.h5 -i excitation/averaged_results.h5 -b 316 -o excitation/band_population --pdf --csv

Default output: PNG + JSON. PDF and CSV are opt-in.
-b only checks the configured initial band; it does not change the data.
""")
    parser.add_argument('-H', '--hamil', required=True, help='HAMIL actually used by this run')
    parser.add_argument('-i', '--input', default='averaged_results.h5', help='Population result HDF5 file')
    parser.add_argument('--dataset', default='sh_pops')
    parser.add_argument('--time-axis', type=int, choices=[0, 1])
    parser.add_argument('--time-dataset')
    parser.add_argument('--dt', type=float, help='Verified output spacing in fs, only if no time dataset exists')
    parser.add_argument('-b', '--init-band', '--expected-initial-band', dest='expected_initial_band', type=int, metavar='BAND', help='Expected initial band for checking, e.g. 316')
    parser.add_argument('-c', '--check-only', action='store_true', help='Print checks without writing files')
    parser.add_argument('-o', '--output', default='band_population', help='Output stem (default: band_population)')
    parser.add_argument('--dpi', type=int, default=600)
    parser.add_argument('--pdf', action='store_true', help='Also save PDF (default: off)')
    parser.add_argument('--csv', action='store_true', help='Also save CSV (default: off)')
    args = parser.parse_args()
    try:
        if args.dt is not None and (not np.isfinite(args.dt) or args.dt <= 0):
            raise ValueError('--dt must be finite and positive.')
        if args.dpi <= 0:
            raise ValueError('--dpi must be positive.')
        hp, rp = Path(args.hamil).resolve(), Path(args.input).resolve()
        print('HAMIL:', hp, flush=True)
        print('Result:', rp, flush=True)
        meta = {}
        with h5py.File(hp, 'r') as hf:
            bands = read_basis(find_dataset(hf, 'basis_list')[...])
            eig = find_dataset(hf, 'eig_t')
            print('HAMIL basis_list:', bands.tolist(), flush=True)
            print('HAMIL eig_t shape:', eig.shape, flush=True)
            if eig.ndim != 2 or eig.shape[1] != len(bands):
                raise ValueError('HAMIL eig_t does not match its own basis_list.')
            for key in ['potim', 'nsw', 'scissor', 'reorder', 'temperature']:
                ds = find_dataset(hf, key, optional=True)
                if ds is not None and ds.size <= 20:
                    meta[key] = np.asarray(ds[...]).tolist()
                    print('HAMIL {}: {}'.format(key, meta[key]))
            ef = find_dataset(hf, 'efield', optional=True)
            if ef is not None:
                if ef.size > 100000:
                    print('HAMIL efield is large; encoding/content needs separate inspection:', ef.shape)
                else:
                    meta['efield'] = decode_efield(ef[...])
                    print('HAMIL embedded efield:\n' + meta['efield'])
        with h5py.File(rp, 'r') as rf:
            ds = find_dataset(rf, args.dataset)
            raw_shape, data_path = list(ds.shape), ds.name
            print('Result {} raw shape: {}'.format(data_path, raw_shape), flush=True)
            rb = find_dataset(rf, 'basis_list', optional=True)
            if rb is not None:
                if not np.array_equal(read_basis(rb[...]), bands):
                    raise ValueError('Result basis_list differs from HAMIL, including column order.')
                print('Result basis_list agrees with HAMIL.')
            else:
                print('No result basis_list: column order is mapped using the selected HAMIL; run provenance remains to be checked.')
            pops, axis = orient_population(ds[...], len(bands), args.time_axis)
            times, time_source = read_time(rf, pops.shape[0], args.time_dataset, args.dt)
        print('Oriented populations: {} samples x {} bands; raw time axis={}'.format(*pops.shape, axis))
        print('Time:', time_source)
        if times is not None:
            print('First/last saved time (fs):', times[0], times[-1])
        print('column0  band      first_saved     last_saved        maximum')
        for col, band in enumerate(bands):
            print('{:7d} {:5d} {:16.9f} {:14.9f} {:14.9f}'.format(col, band, pops[0, col], pops[-1, col], pops[:, col].max()))
        initial_band = int(bands[np.argmax(pops[0])])
        print('Largest first-sample population: band', initial_band)
        if args.expected_initial_band is not None:
            if args.expected_initial_band not in bands:
                raise ValueError('Expected initial band is absent from HAMIL.')
            print('Configured initial band supplied for comparison:', args.expected_initial_band)
            idx = int(np.flatnonzero(bands == args.expected_initial_band)[0])
            if not np.isclose(pops[0, idx], 1., atol=1e-6, rtol=0):
                print('WARNING: first saved population is not a pure configured initial state. Check saved time and run provenance.')
        sums = pops.sum(axis=1)
        norm_error = float(np.max(np.abs(sums - 1)))
        print('Population sum min/max: {:.12f} / {:.12f}'.format(sums.min(), sums.max()))
        print('Maximum normalization error over all saved samples: {:.6g}'.format(norm_error))
        if norm_error > 1e-6 or pops.min() < -1e-8 or pops.max() > 1 + 1e-8:
            print('WARNING: normalization or population bounds need review; data have not been clipped or renormalized.')
        if args.check_only:
            print('CHECK ONLY: no files written; dimensional agreement does not prove common run provenance.')
            return 0
        if times is None:
            raise ValueError('No time coordinate: use --check-only or supply a verified --dt. No plot written.')
        stem = Path(args.output)
        stem.parent.mkdir(parents=True, exist_ok=True)
        plot_population(times, bands, pops, stem, args.dpi, save_pdf=args.pdf)
        if args.csv:
            np.savetxt(str(stem) + '.csv', np.column_stack([times, pops]), delimiter=',',
                       header='time_fs,' + ','.join('band_{}'.format(b) for b in bands), comments='', fmt='%.12g')
        report = dict(hamil=str(hp), result=str(rp), basis_list=bands.tolist(),
                      dataset=data_path, raw_shape=raw_shape, raw_time_axis=axis,
                      time_source=time_source, first_saved=pops[0].tolist(),
                      last_saved=pops[-1].tolist(), normalization_error=norm_error,
                      hamil_metadata=meta, expected_initial_band=args.expected_initial_band,
                      provenance_note='Selected by caller; shape agreement does not establish a common run.')
        Path(str(stem) + '.json').write_text(json.dumps(report, indent=2, ensure_ascii=False), encoding='utf-8')
        saved = [str(stem) + '.png', str(stem) + '.json']
        if args.pdf:
            saved.append(str(stem) + '.pdf')
        if args.csv:
            saved.append(str(stem) + '.csv')
        print('Saved:', ', '.join(saved))
        return 0
    except (OSError, KeyError, ValueError, ImportError) as exc:
        print('ERROR:', exc, file=sys.stderr)
        return 1


if __name__ == '__main__':
    raise SystemExit(main())
