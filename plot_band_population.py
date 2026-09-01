#!/usr/bin/env python3
"""Plot band populations from a surfhop averaged_results.h5 file.

Put this script in the ``excitation`` directory and run:

    python plot_band_population.py

The default labels are bands 316--325.  For another basis, use for example:

    python plot_band_population.py --bands 300-340
    python plot_band_population.py --first-band 310

Outputs:
    band_population.png   high-resolution raster figure
    band_population.pdf   vector figure for slides/papers
    band_population.csv   plotted numerical data
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import h5py
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Plot NAMD-LMI surface-hopping band populations."
    )
    parser.add_argument(
        "--input",
        default="averaged_results.h5",
        help="input HDF5 file (default: averaged_results.h5)",
    )
    parser.add_argument(
        "--dataset",
        default="sh_pops",
        help="population dataset name (default: sh_pops)",
    )
    parser.add_argument(
        "--bands",
        default=None,
        help="inclusive band range/list, e.g. 316-325 or 316,317,320",
    )
    parser.add_argument(
        "--first-band",
        type=int,
        default=316,
        help="label of the first column when --bands is omitted (default: 316)",
    )
    parser.add_argument(
        "--highlight",
        default="316,317",
        help="bands drawn with thicker lines (default: 316,317; use 'none' to disable)",
    )
    parser.add_argument(
        "--dt",
        type=float,
        default=1.0,
        help="time step in fs if no time dataset is present (default: 1.0)",
    )
    parser.add_argument(
        "--output",
        default="band_population",
        help="output filename stem (default: band_population)",
    )
    parser.add_argument(
        "--dpi",
        type=int,
        default=600,
        help="PNG resolution (default: 600)",
    )
    parser.add_argument(
        "--no-pdf",
        action="store_true",
        help="do not create the vector PDF copy",
    )
    return parser.parse_args()


def dataset_paths(handle: h5py.File) -> list[str]:
    paths: list[str] = []

    def visitor(name: str, obj: h5py.Dataset) -> None:
        if isinstance(obj, h5py.Dataset):
            paths.append(name)

    handle.visititems(visitor)
    return paths


def find_dataset(handle: h5py.File, requested: str) -> h5py.Dataset:
    requested = requested.lstrip("/")
    if requested in handle and isinstance(handle[requested], h5py.Dataset):
        return handle[requested]

    matches = [p for p in dataset_paths(handle) if p.split("/")[-1] == requested]
    if len(matches) == 1:
        return handle[matches[0]]
    if not matches:
        available = ", ".join("/" + p for p in dataset_paths(handle)) or "(none)"
        raise KeyError(
            f"dataset '{requested}' was not found. Available datasets: {available}"
        )
    raise KeyError(
        f"dataset name '{requested}' is ambiguous: "
        + ", ".join("/" + p for p in matches)
    )


def parse_integer_list(spec: str) -> list[int]:
    values: list[int] = []
    for token in spec.replace(" ", "").split(","):
        if not token:
            continue
        if "-" in token[1:]:
            split_at = token[1:].index("-") + 1
            start = int(token[:split_at])
            stop = int(token[split_at + 1 :])
            step = 1 if stop >= start else -1
            values.extend(range(start, stop + step, step))
        elif ":" in token:
            start_text, stop_text = token.split(":", maxsplit=1)
            start, stop = int(start_text), int(stop_text)
            step = 1 if stop >= start else -1
            values.extend(range(start, stop + step, step))
        else:
            values.append(int(token))
    if not values:
        raise ValueError("the band list is empty")
    return values


def orient_population(population: np.ndarray, expected_bands: int | None) -> np.ndarray:
    population = np.squeeze(population)
    if population.ndim != 2:
        raise ValueError(
            f"population dataset must become 2-D after squeeze; got shape {population.shape}"
        )

    if expected_bands is not None:
        row_match = population.shape[0] == expected_bands
        col_match = population.shape[1] == expected_bands
        if col_match and not row_match:
            return population
        if row_match and not col_match:
            return population.T
        if not row_match and not col_match:
            raise ValueError(
                f"--bands gives {expected_bands} labels, but dataset shape is "
                f"{population.shape}; neither axis has that length"
            )

    # NAMD trajectories normally contain many more time points than bands.
    if population.shape[0] < population.shape[1]:
        return population.T
    return population


def read_time(handle: h5py.File, ntime: int, dt: float) -> tuple[np.ndarray, str]:
    candidates = ("time", "times", "time_fs", "t")
    paths = dataset_paths(handle)
    for candidate in candidates:
        matches = [p for p in paths if p.split("/")[-1].lower() == candidate]
        for path in matches:
            values = np.asarray(handle[path][...]).squeeze()
            if values.ndim == 1 and values.size == ntime:
                return values.astype(float), "/" + path
    return np.arange(ntime, dtype=float) * dt, f"generated with dt={dt:g} fs"


def finite_real_array(values: np.ndarray) -> np.ndarray:
    if np.iscomplexobj(values):
        max_imag = float(np.nanmax(np.abs(values.imag)))
        if max_imag > 1.0e-10:
            print(
                f"Warning: population has a non-negligible imaginary part "
                f"(max |Im|={max_imag:.3e}); plotting the real part.",
                file=sys.stderr,
            )
        values = values.real
    values = np.asarray(values, dtype=float)
    if not np.all(np.isfinite(values)):
        raise ValueError("population dataset contains NaN or infinite values")
    return values


def save_csv(path: Path, time_fs: np.ndarray, bands: list[int], population: np.ndarray) -> None:
    table = np.column_stack((time_fs, population))
    header = "time_fs," + ",".join(f"band_{band}" for band in bands)
    np.savetxt(path, table, delimiter=",", header=header, comments="", fmt="%.10g")


def plot_population(
    time_fs: np.ndarray,
    bands: list[int],
    population: np.ndarray,
    highlighted: set[int],
    png_path: Path,
    pdf_path: Path | None,
    dpi: int,
) -> None:
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 11,
            "axes.linewidth": 1.0,
            "xtick.direction": "in",
            "ytick.direction": "in",
            "xtick.top": True,
            "ytick.right": True,
        }
    )
    fig, ax = plt.subplots(figsize=(8.2, 5.0))
    colors = plt.get_cmap("tab10")(np.linspace(0.0, 1.0, max(10, len(bands))))

    for index, band in enumerate(bands):
        is_main = band in highlighted
        ax.plot(
            time_fs,
            population[:, index],
            color=colors[index % len(colors)],
            linewidth=2.4 if is_main else 1.25,
            alpha=1.0 if is_main else 0.78,
            label=f"Band {band}",
            zorder=3 if is_main else 2,
        )

    ymin = min(0.0, float(np.min(population)))
    ymax = max(1.0, float(np.max(population)))
    margin = max(0.02, 0.03 * (ymax - ymin if ymax > ymin else 1.0))
    ax.set_xlim(float(time_fs[0]), float(time_fs[-1]))
    ax.set_ylim(ymin - margin, ymax + margin)
    ax.set_xlabel("Time (fs)")
    ax.set_ylabel("Population")
    ax.grid(True, linestyle="--", linewidth=0.55, alpha=0.32)
    ax.legend(
        loc="center left",
        bbox_to_anchor=(1.01, 0.5),
        frameon=False,
        ncol=1 if len(bands) <= 12 else 2,
        fontsize=9,
    )
    fig.tight_layout()
    fig.savefig(png_path, dpi=dpi, bbox_inches="tight")
    if pdf_path is not None:
        fig.savefig(pdf_path, bbox_inches="tight")
    plt.close(fig)


def main() -> int:
    args = parse_args()
    input_path = Path(args.input)
    if not input_path.is_file():
        print(
            f"Error: '{input_path}' was not found. Run this script inside the "
            "surfhop-generated excitation directory, or pass --input FILE.",
            file=sys.stderr,
        )
        return 1
    if args.dt <= 0:
        print("Error: --dt must be positive.", file=sys.stderr)
        return 1

    try:
        explicit_bands = parse_integer_list(args.bands) if args.bands else None
        with h5py.File(input_path, "r") as handle:
            dataset = find_dataset(handle, args.dataset)
            dataset_name = dataset.name
            population = orient_population(
                finite_real_array(dataset[...]),
                len(explicit_bands) if explicit_bands is not None else None,
            )
            time_fs, time_source = read_time(handle, population.shape[0], args.dt)

        if explicit_bands is None:
            bands = list(range(args.first_band, args.first_band + population.shape[1]))
        else:
            bands = explicit_bands
        if len(bands) != population.shape[1]:
            raise ValueError(
                f"there are {population.shape[1]} population columns but "
                f"{len(bands)} band labels"
            )

        if args.highlight.strip().lower() in {"", "none", "off"}:
            highlighted: set[int] = set()
        else:
            highlighted = set(parse_integer_list(args.highlight))

        output_stem = Path(args.output)
        png_path = output_stem.with_suffix(".png")
        pdf_path = None if args.no_pdf else output_stem.with_suffix(".pdf")
        csv_path = output_stem.with_suffix(".csv")

        save_csv(csv_path, time_fs, bands, population)
        plot_population(
            time_fs,
            bands,
            population,
            highlighted,
            png_path,
            pdf_path,
            args.dpi,
        )

    except (OSError, KeyError, ValueError) as exc:
        print(f"Error: {exc}", file=sys.stderr)
        return 1

    print(f"Read {dataset_name}: {population.shape[0]} time points x {population.shape[1]} bands")
    print(f"Time axis: {time_source}")
    print("Population summary:")
    for index, band in enumerate(bands):
        values = population[:, index]
        imax = int(np.argmax(values))
        print(
            f"  band {band}: initial={values[0]:.6f}, final={values[-1]:.6f}, "
            f"max={values[imax]:.6f} at {time_fs[imax]:.1f} fs"
        )

    sums = np.sum(population, axis=1)
    print(
        f"Population sum: initial={sums[0]:.6f}, final={sums[-1]:.6f}, "
        f"range=[{np.min(sums):.6f}, {np.max(sums):.6f}]"
    )
    print(f"Saved: {png_path}")
    if pdf_path is not None:
        print(f"Saved: {pdf_path}")
    print(f"Saved: {csv_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
