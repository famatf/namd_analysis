#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Extract selected bands from a continuous NVE OUTCAR and optionally run FFT.

Dependencies: Python >= 3.9, numpy; matplotlib only for PNG/PDF.
One complete selected k-point/spin table per E-fermi block is retained.
This assumes one final eigenvalue block per MD frame; SCF convergence and
omitted MD frames are not verified. An incomplete interior table is rejected.
Band numbers are fixed VASP indices, not dynamically tracked orbital character.
No scissor correction, smoothing or interpolation is applied.
"""
from __future__ import annotations
import argparse
import json
import re
import sys
from dataclasses import dataclass
from pathlib import Path
import numpy as np

FLOAT_TEXT = r"[-+]?(?:(?:\d+(?:\.\d*)?)|(?:\.\d+))(?:[EeDd][-+]?\d+)?"
FERMI_RE = re.compile(rf"E-fermi\s*:\s*({FLOAT_TEXT})")
POTIM_RE = re.compile(rf"\bPOTIM\s*=\s*({FLOAT_TEXT})")
NBANDS_RE = re.compile(r"\bNBANDS\s*=\s*(\d+)")
SPIN_RE = re.compile(r"^\s*spin component\s+(\d+)", re.IGNORECASE)
KPOINT_RE = re.compile(r"^\s*k-point\s+(\d+)\s*:", re.IGNORECASE)
BAND_HEADER_RE = re.compile(
    r"^\s*band\s+No\.\s+band\s+energies(?:\s+occupation)?",
    re.IGNORECASE,
)
BAND_ROW_RE = re.compile(
    rf"^\s*(\d+)\s+({FLOAT_TEXT})(?:\s+({FLOAT_TEXT}))?\s*$"
)


@dataclass
class StepResult:
    """One completed NVE ionic step extracted from the MD OUTCAR."""

    source_table: int
    efermi: float
    energies: dict[int, float]
    occupations: dict[int, float]


@dataclass
class ParseResult:
    """All usable NVE steps and metadata found in one OUTCAR."""

    steps: list[StepResult]
    potim_fs: float | None
    nbands: int | None
    ignored_trailing_table: bool


def as_float(text: str) -> float:
    return float(text.replace("D", "E").replace("d", "e"))


def read_nve_outcar(
    path: Path,
    requested_bands: tuple[int, ...],
    requested_kpoint: int,
    requested_spin: int,
) -> ParseResult:
    """
    Extract one target eigenvalue table for each completed MD ionic step.

    VASP prints an E-fermi line followed by the k-point band tables after the
    electronic minimization of an ionic step.  The parser associates a target
    table with the most recent E-fermi line.  A complete target table becomes
    one NVE frame.

    If the file is being written while the calculation runs, an incomplete
    final table is ignored.  An incomplete target table in the middle of the
    trajectory is treated as a hard error because silently skipping it would
    make the time series discontinuous.
    """
    if not path.is_file():
        raise FileNotFoundError(f"Cannot find NVE OUTCAR: {path}")

    requested_set = set(requested_bands)
    steps: list[StepResult] = []
    incomplete_target_tables: list[int] = []

    nbands: int | None = None
    potim_fs: float | None = None
    current_spin = 1
    pending_kpoint: int | None = None
    latest_fermi: float | None = None
    fermi_serial = 0

    table_active = False
    table_spin = 1
    table_kpoint = -1
    table_band_numbers: set[int] = set()
    table_energies: dict[int, float] = {}
    table_occupations: dict[int, float] = {}
    target_table_serial = 0

    # A target table normally occurs once after each E-fermi.  Retaining this
    # serial prevents an accidental duplicate table from becoming a fake MD
    # step; the later complete table replaces the earlier one.
    last_record_fermi_serial: int | None = None

    def finish_table() -> None:
        nonlocal table_active
        nonlocal pending_kpoint
        nonlocal table_band_numbers
        nonlocal table_energies
        nonlocal table_occupations
        nonlocal target_table_serial
        nonlocal last_record_fermi_serial

        if not table_active:
            return

        is_target = (
            table_spin == requested_spin
            and table_kpoint == requested_kpoint
        )

        if is_target:
            target_table_serial += 1
            complete_by_count = (
                nbands is None or table_band_numbers == set(range(1, nbands + 1))
            )
            complete_by_selection = requested_set.issubset(table_energies)

            if complete_by_count and complete_by_selection:
                if latest_fermi is None:
                    raise ValueError(
                        "Found a target band table before any E-fermi line."
                    )

                result = StepResult(
                    source_table=target_table_serial,
                    efermi=latest_fermi,
                    energies=dict(table_energies),
                    occupations=dict(table_occupations),
                )

                if (
                    steps
                    and last_record_fermi_serial == fermi_serial
                ):
                    steps[-1] = result
                else:
                    steps.append(result)
                last_record_fermi_serial = fermi_serial
            else:
                incomplete_target_tables.append(target_table_serial)

        table_active = False
        pending_kpoint = None
        table_band_numbers = set()
        table_energies = {}
        table_occupations = {}

    with path.open("r", encoding="utf-8", errors="ignore") as handle:
        for line in handle:
            if table_active:
                row_match = BAND_ROW_RE.match(line)
                if row_match:
                    band_number = int(row_match.group(1))
                    table_band_numbers.add(band_number)
                    if band_number in requested_set:
                        table_energies[band_number] = as_float(
                            row_match.group(2)
                        )
                        occupation_text = row_match.group(3)
                        table_occupations[band_number] = (
                            as_float(occupation_text)
                            if occupation_text is not None
                            else float("nan")
                        )
                    continue

                # Empty lines before the first band row are harmless.  Once at
                # least one band row has appeared, the first non-row closes the
                # table.  The same line is then processed for other markers.
                if table_band_numbers or line.strip():
                    finish_table()

            nbands_match = NBANDS_RE.search(line)
            if nbands_match:
                nbands = int(nbands_match.group(1))

            if potim_fs is None:
                potim_match = POTIM_RE.search(line)
                if potim_match:
                    potim_fs = as_float(potim_match.group(1))

            fermi_match = FERMI_RE.search(line)
            if fermi_match:
                latest_fermi = as_float(fermi_match.group(1))
                fermi_serial += 1
                current_spin = 1

            spin_match = SPIN_RE.match(line)
            if spin_match:
                current_spin = int(spin_match.group(1))
                continue

            kpoint_match = KPOINT_RE.match(line)
            if kpoint_match:
                pending_kpoint = int(kpoint_match.group(1))
                continue

            if BAND_HEADER_RE.match(line) and pending_kpoint is not None:
                table_active = True
                table_spin = current_spin
                table_kpoint = pending_kpoint
                table_band_numbers = set()
                table_energies = {}
                table_occupations = {}

    finish_table()

    if not steps:
        raise ValueError(
            "No complete target band tables were found. Check that this is "
            "the NVE OUTCAR, that VASP printed eigenvalues, and that --kpoint, "
            "--spin, --bands, and NBANDS are correct."
        )

    ignored_trailing_table = False
    if incomplete_target_tables:
        last_incomplete = incomplete_target_tables[-1]
        incomplete_before_last_complete = [
            number
            for number in incomplete_target_tables
            if number < steps[-1].source_table
        ]
        if incomplete_before_last_complete:
            preview = ", ".join(
                str(number) for number in incomplete_before_last_complete[:10]
            )
            raise ValueError(
                "Incomplete target band table(s) occurred inside the NVE "
                f"trajectory: {preview}. Refusing to skip an interior frame."
            )
        if last_incomplete > steps[-1].source_table:
            ignored_trailing_table = True

    return ParseResult(
        steps=steps,
        potim_fs=potim_fs,
        nbands=nbands,
        ignored_trailing_table=ignored_trailing_table,
    )


def spectrum(values, dt, window='hann', detrend='mean'):
    """One-sided FFT amplitude in eV; never independently normalize curves."""
    y = np.asarray(values, dtype=float)
    if y.ndim == 1:
        y = y[:, None]
    n = len(y)
    if n < 4:
        raise ValueError('FFT requires at least four selected frames.')
    y = y - y.mean(axis=0)
    if detrend == 'linear':
        x = np.arange(n, dtype=float) - (n - 1) / 2
        y = y - x[:, None] * (x @ y / (x @ x))[None, :]
    w = np.hanning(n) if window == 'hann' else np.ones(n)
    amp = np.abs(np.fft.rfft(y * w[:, None], axis=0)) / w.sum()
    amp[1:] *= 2
    if n % 2 == 0:
        amp[-1] /= 2
    freq = np.fft.rfftfreq(n, d=dt)
    # Speed of light is exactly 299792458 m/s; input time is fs.
    wn = freq * 1e15 / 2.99792458e10
    return freq * 1000, wn, amp


def arguments():
    p = argparse.ArgumentParser(description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='''Examples (run in the NVE directory):
  Band energies only; default output is PNG:
    python nve_band_analysis.py -b 313 316 317 321
  Individual-band FFT:
    python nve_band_analysis.py -b 313 316 317 321 --fft
  Gap E317-E316 and its FFT, without individual-band FFT:
    python nve_band_analysis.py -b 316 317 -g 316 317 --gap-fft
  Both FFT types; export every format:
    python nve_band_analysis.py -i OUTCAR -b 313 316 317 321 -f -g 316 317 --gap-fft --save png pdf csv json -o nve
  Gap trace without FFT; selected frames; CSV and JSON only:
    python nve_band_analysis.py -b 316 317 -g 316 317 --start 501 --end 2500 --save csv json

-g A B always means E_B - E_A, not an absolute difference.
Gap bands are read automatically even if absent from -b.
--fft and --gap-fft are independent, and both default to off.
--save replaces the format list; default: png. Other formats are opt-in.
FFT uses equally spaced frames, mean removal, and a Hann window by default.
Frequency units: cm^-1 (plot), with THz also included in CSV.
''')
    p.add_argument('-i', '--outcar', default='OUTCAR')
    p.add_argument('-b', '--bands', nargs='+', type=int, required=True, help='Band IDs, separated by spaces')
    p.add_argument('-f', '--fft', action='store_true', help='FFT of every band selected with -b')
    p.add_argument('-g', '--gap', nargs=2, type=int, metavar=('A','B'), help='Analyze E_B - E_A')
    p.add_argument('--gap-fft', action='store_true', help='FFT of the gap specified by -g')
    p.add_argument('--save', '--formats', nargs='+', choices=['png','pdf','csv','json'], default=['png'], help='Output formats (default: png)')
    p.add_argument('-o', '--prefix', default='nve_bands', help='Output prefix, optionally with directory')
    p.add_argument('--start', type=int, default=1, help='First retained frame, one-based inclusive')
    p.add_argument('--end', type=int, help='Last retained frame, inclusive')
    p.add_argument('--dt', type=float, help='Saved-frame spacing in fs; otherwise OUTCAR POTIM')
    p.add_argument('-k', '--kpoint', type=int, default=1)
    p.add_argument('-s', '--spin', type=int, default=1)
    p.add_argument('--reference', choices=['mean-fermi','raw'], default='mean-fermi', help='Shared fixed reference for energy traces')
    p.add_argument('--window', choices=['hann','none'], default='hann')
    p.add_argument('--detrend', choices=['mean','linear'], default='mean')
    p.add_argument('--fmax', type=float, help='FFT plot upper limit in cm^-1; exported data remain complete')
    p.add_argument('--dpi', type=int, default=400)
    return p.parse_args()


def main():
    args = arguments()
    bands = tuple(args.bands)
    if len(set(bands)) != len(bands) or any(b < 1 for b in bands):
        raise ValueError('Band IDs must be unique positive integers.')
    if args.gap and (min(args.gap) < 1 or args.gap[0] == args.gap[1]):
        raise ValueError('Gap requires two distinct positive band IDs.')
    if args.gap_fft and not args.gap:
        raise ValueError('--gap-fft requires -g A B.')
    if args.start < 1 or (args.end is not None and args.end < args.start):
        raise ValueError('Require 1 <= start <= end.')
    if min(args.spin,args.kpoint,args.dpi) < 1:
        raise ValueError('Spin, k-point and DPI must be positive.')
    if args.fmax is not None and (not np.isfinite(args.fmax) or args.fmax <= 0):
        raise ValueError('--fmax must be finite and positive.')
    read_bands = tuple(dict.fromkeys([*bands, *(args.gap or [])]))
    source = Path(args.outcar).resolve()
    parsed = read_nve_outcar(source, read_bands, args.kpoint, args.spin)
    end = args.end if args.end is not None else len(parsed.steps)
    if end > len(parsed.steps):
        raise ValueError('Requested end={} but only {} complete retained frames exist.'.format(end,len(parsed.steps)))
    steps = parsed.steps[args.start-1:end]
    if len(steps) < 2:
        raise ValueError('Select at least two complete frames.')
    dt = args.dt if args.dt is not None else parsed.potim_fs
    if dt is None or not np.isfinite(dt) or dt <= 0:
        raise ValueError('Need positive finite OUTCAR POTIM or explicit --dt.')
    time = np.arange(len(steps),dtype=float)*dt
    raw = np.array([[s.energies[b] for b in read_bands] for s in steps])
    fermi = np.array([s.efermi for s in steps])
    if not np.all(np.isfinite(raw)) or not np.all(np.isfinite(fermi)):
        raise ValueError('Non-finite energies or Fermi values.')
    reference = float(fermi.mean()) if args.reference == 'mean-fermi' else 0.
    idx = {b:i for i,b in enumerate(read_bands)}
    plotted = raw[:,[idx[b] for b in bands]]-reference
    gap = raw[:,idx[args.gap[1]]]-raw[:,idx[args.gap[0]]] if args.gap else None
    specs = {}
    if args.fft:
        specs['band_fft'] = spectrum(plotted,dt,args.window,args.detrend)
    if args.gap_fft:
        specs['gap_fft'] = spectrum(gap,dt,args.window,args.detrend)
    prefix = Path(args.prefix)
    prefix.parent.mkdir(parents=True,exist_ok=True)
    formats = set(args.save)
    saved = []
    def filename(tag, ext):
        return Path(str(prefix)+'_'+tag+'.'+ext)
    def csv(tag, cols, names):
        path=filename(tag,'csv')
        np.savetxt(path,np.column_stack(cols),delimiter=',',header=','.join(names),comments='',fmt='%.12g')
        saved.append(str(path))
    def chart(tag,x,y,labels,xlabel,ylabel):
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        y=np.asarray(y)
        if y.ndim==1: y=y[:,None]
        fig,ax=plt.subplots(figsize=(10,6))
        palette=plt.get_cmap('tab10').colors if y.shape[1]<=10 else plt.get_cmap('tab20').colors
        for c,label in enumerate(labels):
            color=palette[c] if c<len(palette) else plt.get_cmap('hsv')(c/y.shape[1])
            ax.plot(x,y[:,c],label=label,color=color,lw=1.3)
        ax.set_xlabel(xlabel); ax.set_ylabel(ylabel)
        ax.set_xlim(float(x[0]),float(x[-1]))
        if tag.endswith('fft') and args.fmax is not None:
            ax.set_xlim(0,min(args.fmax,float(x[-1])))
        ax.grid(alpha=.25)
        ax.legend(loc='upper center',bbox_to_anchor=(.5,1.15),ncol=min(4,len(labels)),frameon=False)
        fig.tight_layout()
        for ext in ('png','pdf'):
            if ext in formats:
                path=filename(tag,ext); fig.savefig(path,dpi=args.dpi,bbox_inches='tight'); saved.append(str(path))
        plt.close(fig)
    if formats & {'png','pdf'}:
        chart('energies',time,plotted,['Band {}'.format(b) for b in bands],'Time (fs)',
              'Energy - mean Fermi energy (eV)' if args.reference=='mean-fermi' else 'Energy (eV; as stored)')
        if gap is not None:
            chart('gap',time,gap,['E{} - E{}'.format(args.gap[1],args.gap[0])],'Time (fs)','Energy difference (eV)')
        for tag,(_,wn,amp) in specs.items():
            labels=['Band {}'.format(b) for b in bands] if tag=='band_fft' else ['E{} - E{}'.format(args.gap[1],args.gap[0])]
            chart(tag,wn,amp,labels,'Wavenumber (cm$^{-1}$)','One-sided FFT amplitude (eV)')
    if 'csv' in formats:
        cols=[np.arange(args.start,end+1),time,fermi]; names=['retained_frame','time_fs','Efermi_eV']
        for b in read_bands:
            cols.extend([raw[:,idx[b]],raw[:,idx[b]]-reference,[s.occupations[b] for s in steps]])
            names.extend(['E{}_raw_eV'.format(b),'E{}_referenced_eV'.format(b),'occupation_{}'.format(b)])
        csv('energies',cols,names)
        if gap is not None: csv('gap',[time,gap],['time_fs','E{}_minus_E{}_eV'.format(args.gap[1],args.gap[0])])
        for tag,(thz,wn,amp) in specs.items():
            labels=['band_{}_amplitude_eV'.format(b) for b in bands] if tag=='band_fft' else ['gap_amplitude_eV']
            csv(tag,[thz,wn,*amp.T],['frequency_THz','wavenumber_cm-1',*labels])
    fft_meta=None
    if specs:
        fft_meta=dict(window=args.window,detrend=args.detrend,normalization='One-sided amplitude / sum(window); double interior bins only',
                      amplitude_unit='eV',bin_spacing_cm1=1e15/(len(time)*dt*2.99792458e10),
                      nyquist_cm1=1e15/(2*dt*2.99792458e10),zero_padding=False,
                      note='Bin spacing is not peak-location uncertainty; windowing affects spectral resolution. No mode assignment is inferred.')
    report=dict(outcar=str(source),bands=list(bands),read_bands=list(read_bands),gap_bands=args.gap,
                gap_definition='E_B - E_A',spin=args.spin,kpoint=args.kpoint,available_frames=len(parsed.steps),
                start=args.start,end=end,selected_frames=len(time),dt_fs=dt,
                time_source='--dt' if args.dt is not None else 'POTIM',time_origin='First selected retained table is 0 fs',
                span_fs=float(time[-1]),reference=args.reference,subtracted_constant_eV=reference,
                band_fft=args.fft,gap_fft=args.gap_fft,fft=fft_meta,
                gap_statistics_eV=None if gap is None else dict(mean=float(gap.mean()),std=float(gap.std()),min=float(gap.min()),max=float(gap.max()),peak_to_peak=float(np.ptp(gap))),
                ignored_incomplete_tail=parsed.ignored_trailing_table,
                assumption='One final selected eigenvalue table per MD frame; complete tables do not verify SCF convergence or omitted MD steps.',
                band_identity='Fixed band index; no layer/orbital tracking',scissor_applied=False,outputs=saved.copy())
    if 'json' in formats:
        path=filename('analysis','json'); report['outputs'].append(str(path))
        path.write_text(json.dumps(report,indent=2,allow_nan=False),encoding='utf-8'); saved.append(str(path))
    print('Input:',source)
    print('Bands:',list(bands),' Gap (A,B):',args.gap)
    print('Frames {}..{}; dt={} fs; relative time=0..{} fs'.format(args.start,end,dt,time[-1]))
    print('Band FFT:',args.fft,' Gap FFT:',args.gap_fft)
    if fft_meta: print('FFT bin spacing: {:.6g} cm^-1; Nyquist: {:.6g} cm^-1'.format(fft_meta['bin_spacing_cm1'],fft_meta['nyquist_cm1']))
    if gap is not None: print('Gap mean/std: {:.9f} / {:.9f} eV'.format(gap.mean(),gap.std()))
    if parsed.ignored_trailing_table: print('WARNING: incomplete trailing target table ignored.')
    print('Time assumes one saved final band table per MD frame; missing frames and SCF convergence are not verified.')
    for path in saved: print('Saved:',path)


if __name__=='__main__':
    try:
        main()
    except (OSError,ValueError,ImportError) as exc:
        print('ERROR:',exc,file=sys.stderr)
        sys.exit(1)


