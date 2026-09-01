#!/usr/bin/env python3
"""Plot layer- and direction-projected Gamma-point phonon modes.

Input is a Phonopy qpoints.yaml containing eigenvectors at q = (0, 0, 0).
If qpoints.yaml omits atom symbols, as some Phonopy versions do, symbols are
read from the companion phonopy.yaml. The script is tailored to a BP/MoS2
structure (P atoms form BP; Mo and S atoms form MoS2), but it does not assume
a particular atom ordering.

Outputs, using the chosen prefix:
  *_projection.png  publication-resolution raster figure
  *_projection.pdf  vector figure
  *_projection.csv  numerical frequencies and projection weights
"""

from __future__ import annotations

import argparse
import csv
import sys
from collections import Counter
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, Normalize
import numpy as np
import yaml


THZ_TO_CM1 = 33.35640951981521


def parse_mode_list(text: str) -> list[int]:
    """Parse comma-separated one-based mode numbers."""
    if not text.strip():
        return []
    try:
        modes = [int(item.strip()) for item in text.split(",") if item.strip()]
    except ValueError as exc:
        raise argparse.ArgumentTypeError(
            "mode numbers must be comma-separated integers"
        ) from exc
    if any(mode < 1 for mode in modes):
        raise argparse.ArgumentTypeError("mode numbers must be positive")
    return modes


def _atom_symbols_from_list(points, natom: int) -> np.ndarray | None:
    """Return symbols from a Phonopy points/atoms list when compatible."""
    if not isinstance(points, list) or len(points) != natom:
        return None
    if not all(isinstance(point, dict) and "symbol" in point for point in points):
        return None
    return np.asarray([point["symbol"] for point in points], dtype=str)


def _find_symbol_list_recursively(node, natom: int) -> np.ndarray | None:
    """Fallback search for a compatible atom list in a Phonopy YAML tree."""
    if isinstance(node, dict):
        for key in ("points", "atoms"):
            symbols = _atom_symbols_from_list(node.get(key), natom)
            if symbols is not None:
                return symbols
        for value in node.values():
            symbols = _find_symbol_list_recursively(value, natom)
            if symbols is not None:
                return symbols
    elif isinstance(node, list):
        for value in node:
            symbols = _find_symbol_list_recursively(value, natom)
            if symbols is not None:
                return symbols
    return None


def load_atom_symbols(qpoint_data, natom: int, structure_yaml: Path) -> tuple[np.ndarray, str]:
    """Load atom symbols from qpoints.yaml or a companion phonopy.yaml."""
    for key in ("points", "atoms"):
        symbols = _atom_symbols_from_list(qpoint_data.get(key), natom)
        if symbols is not None:
            return symbols, f"qpoints.yaml:{key}"

    if not structure_yaml.is_file():
        raise ValueError(
            "qpoints.yaml omits atom symbols and the companion structure file "
            f"was not found: {structure_yaml}. Put phonopy.yaml beside "
            "qpoints.yaml or pass --structure-yaml PATH."
        )

    with structure_yaml.open("r", encoding="utf-8") as stream:
        structure_data = yaml.safe_load(stream)

    # Phonopy q-point eigenvectors refer to primitive-cell atoms. Prefer that
    # list, then accept unit/supercell lists only when their length is natom.
    for container_name in ("primitive_cell", "unit_cell", "supercell"):
        container = structure_data.get(container_name)
        if not isinstance(container, dict):
            continue
        for key in ("points", "atoms"):
            symbols = _atom_symbols_from_list(container.get(key), natom)
            if symbols is not None:
                return symbols, f"{structure_yaml.name}:{container_name}/{key}"

    symbols = _find_symbol_list_recursively(structure_data, natom)
    if symbols is not None:
        return symbols, f"{structure_yaml.name}:atom list"

    raise ValueError(
        f"could not find a {natom}-atom symbol list in {structure_yaml}"
    )


def load_gamma_modes(path: Path, structure_yaml: Path):
    """Read symbols, frequencies, and complex eigenvectors at Gamma."""
    with path.open("r", encoding="utf-8") as stream:
        data = yaml.safe_load(stream)

    natom = int(data["natom"])
    symbols, symbol_source = load_atom_symbols(data, natom, structure_yaml)

    unexpected = sorted(set(symbols) - {"P", "Mo", "S"})
    if unexpected:
        raise ValueError(
            "this BP/MoS2 projection script found unexpected elements: "
            + ", ".join(unexpected)
        )

    phonons = data.get("phonon", [])
    if not phonons:
        raise ValueError("qpoints.yaml contains no phonon q-points")

    q_vectors = np.asarray(
        [entry["q-position"] for entry in phonons], dtype=float
    )
    gamma_index = int(np.argmin(np.linalg.norm(q_vectors, axis=1)))
    gamma_q = q_vectors[gamma_index]
    if np.linalg.norm(gamma_q) > 1.0e-7:
        raise ValueError(
            f"no Gamma point found; closest q-point is {gamma_q.tolist()}"
        )

    bands = phonons[gamma_index].get("band", [])
    expected_modes = 3 * natom
    if len(bands) != expected_modes:
        raise ValueError(
            f"Gamma contains {len(bands)} modes; expected 3N={expected_modes}"
        )

    frequencies = np.empty(expected_modes, dtype=float)
    eigenvectors = np.empty((expected_modes, natom, 3), dtype=complex)

    for mode_index, band in enumerate(bands):
        frequencies[mode_index] = float(band["frequency"])
        raw = np.asarray(band["eigenvector"], dtype=float)
        if raw.shape != (natom, 3, 2):
            raise ValueError(
                "unexpected eigenvector shape for mode "
                f"{mode_index + 1}: {raw.shape}; expected {(natom, 3, 2)}"
            )
        eigenvectors[mode_index] = raw[..., 0] + 1j * raw[..., 1]

    return symbols, frequencies, eigenvectors, gamma_q, symbol_source


def calculate_projections(symbols: np.ndarray, eigenvectors: np.ndarray):
    """Calculate normalized |e|^2 layer and Cartesian projections."""
    power = np.abs(eigenvectors) ** 2
    totals = power.sum(axis=(1, 2))
    if np.any(totals <= 0):
        raise ValueError("one or more eigenvectors have zero norm")
    power /= totals[:, None, None]

    bp_mask = symbols == "P"
    mos2_mask = np.isin(symbols, ["Mo", "S"])

    bp = power[:, bp_mask, :].sum(axis=(1, 2))
    mos2 = power[:, mos2_mask, :].sum(axis=(1, 2))
    xyz = power.sum(axis=1)
    x, y, z = xyz[:, 0], xyz[:, 1], xyz[:, 2]
    in_plane = x + y

    if not np.allclose(bp + mos2, 1.0, atol=1.0e-8):
        raise ValueError("BP and MoS2 projection weights do not sum to one")
    if not np.allclose(in_plane + z, 1.0, atol=1.0e-8):
        raise ValueError("in-plane and out-of-plane weights do not sum to one")

    return bp, mos2, x, y, z, in_plane


def write_csv(
    path: Path,
    frequencies_cm1: np.ndarray,
    frequencies_thz: np.ndarray,
    bp: np.ndarray,
    mos2: np.ndarray,
    x: np.ndarray,
    y: np.ndarray,
    z: np.ndarray,
    in_plane: np.ndarray,
) -> None:
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        writer.writerow(
            [
                "mode",
                "frequency_cm-1",
                "frequency_THz",
                "BP_weight",
                "MoS2_weight",
                "x_weight",
                "y_weight",
                "z_weight",
                "in_plane_weight",
            ]
        )
        for row in zip(
            range(1, len(frequencies_cm1) + 1),
            frequencies_cm1,
            frequencies_thz,
            bp,
            mos2,
            x,
            y,
            z,
            in_plane,
        ):
            writer.writerow([row[0], *[f"{value:.10f}" for value in row[1:]]])


def style_axis(ax: plt.Axes, mode_count: int) -> None:
    ax.axhline(0.0, color="#444444", linewidth=0.8, linestyle="--", zorder=1)
    ax.set_xlim(0.5, mode_count + 0.5)
    ax.set_xlabel("Mode index at Gamma")
    ax.set_ylabel(r"Frequency (cm$^{-1}$)")
    ax.grid(axis="y", color="#D9D9D9", linewidth=0.6, alpha=0.7)
    ax.set_axisbelow(True)


def zoom_mode_axis_to_frequency_window(
    ax: plt.Axes,
    frequencies: np.ndarray,
    frequency_min: float,
    frequency_max: float,
) -> None:
    """Limit mode-index x range to modes visible in a frequency zoom."""
    visible_modes = (
        np.flatnonzero(
            (frequencies >= frequency_min) & (frequencies <= frequency_max)
        )
        + 1
    )
    if len(visible_modes):
        ax.set_xlim(
            max(0.5, visible_modes.min() - 1.5),
            min(len(frequencies) + 0.5, visible_modes.max() + 1.5),
        )


def mark_selected_modes(
    ax: plt.Axes,
    modes: np.ndarray,
    frequencies: np.ndarray,
    selected_modes: list[int],
) -> None:
    x_min, x_max = ax.get_xlim()
    y_min, y_max = ax.get_ylim()
    for label_index, mode in enumerate(selected_modes):
        if mode > len(frequencies):
            continue
        frequency = frequencies[mode - 1]
        if mode < min(x_min, x_max) or mode > max(x_min, x_max):
            continue
        if frequency < min(y_min, y_max) or frequency > max(y_min, y_max):
            continue
        ax.scatter(
            [mode],
            [frequency],
            s=78,
            facecolors="none",
            edgecolors="#111111",
            linewidths=1.2,
            zorder=5,
        )
        vertical_offset = 8 if label_index % 2 == 0 else -14
        ax.annotate(
            f"{mode:03d}\n{frequency:.2f}",
            xy=(mode, frequency),
            xytext=(5, vertical_offset),
            textcoords="offset points",
            fontsize=7.5,
            color="#111111",
            ha="left",
            va="bottom" if vertical_offset > 0 else "top",
            zorder=6,
        )


def plot_projection(
    output_prefix: Path,
    frequencies: np.ndarray,
    bp: np.ndarray,
    z: np.ndarray,
    selected_modes: list[int],
    low_max: float,
    high_min: float,
    high_max: float,
    dpi: int,
) -> tuple[Path, Path]:
    modes = np.arange(1, len(frequencies) + 1)
    norm = Normalize(vmin=0.0, vmax=1.0)

    layer_cmap = LinearSegmentedColormap.from_list(
        "MoS2_to_BP", ["#D55E00", "#ECECEC", "#0072B2"]
    )
    direction_cmap = LinearSegmentedColormap.from_list(
        "in_plane_to_out_of_plane", ["#009E73", "#ECECEC", "#CC79A7"]
    )

    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 10,
            "axes.titlesize": 11,
            "axes.labelsize": 10,
            "xtick.labelsize": 8.5,
            "ytick.labelsize": 8.5,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )

    fig, axes = plt.subplots(
        2,
        2,
        figsize=(13.2, 9.4),
        constrained_layout=True,
        gridspec_kw={"height_ratios": [1.08, 1.0]},
    )
    ax_layer, ax_direction, ax_low, ax_high = axes.ravel()

    scatter_options = {
        "s": 23,
        "linewidths": 0.15,
        "edgecolors": "#333333",
        "alpha": 0.96,
        "zorder": 3,
    }

    layer_scatter = ax_layer.scatter(
        modes,
        frequencies,
        c=bp,
        cmap=layer_cmap,
        norm=norm,
        **scatter_options,
    )
    ax_layer.set_title("(a) Layer projection")
    style_axis(ax_layer, len(frequencies))
    layer_bar = fig.colorbar(layer_scatter, ax=ax_layer, pad=0.015)
    layer_bar.set_label(r"BP weight  ($\sum |e|^2$)")
    layer_bar.set_ticks([0.0, 0.5, 1.0])
    layer_bar.set_ticklabels([r"MoS$_2$", "mixed", "BP"])

    direction_scatter = ax_direction.scatter(
        modes,
        frequencies,
        c=z,
        cmap=direction_cmap,
        norm=norm,
        **scatter_options,
    )
    ax_direction.set_title("(b) Polarization projection")
    style_axis(ax_direction, len(frequencies))
    direction_bar = fig.colorbar(direction_scatter, ax=ax_direction, pad=0.015)
    direction_bar.set_label(r"Out-of-plane weight  ($\sum |e_z|^2$)")
    direction_bar.set_ticks([0.0, 0.5, 1.0])
    direction_bar.set_ticklabels(["in-plane", "mixed", "out-of-plane"])

    ax_low.scatter(
        modes,
        frequencies,
        c=bp,
        cmap=layer_cmap,
        norm=norm,
        **scatter_options,
    )
    ax_low.set_title(f"(c) Low-frequency zoom (<= {low_max:g} cm$^{{-1}}$)")
    style_axis(ax_low, len(frequencies))
    low_floor = min(-10.0, 5.0 * np.floor(frequencies.min() / 5.0) - 2.0)
    ax_low.set_ylim(low_floor, low_max)
    zoom_mode_axis_to_frequency_window(ax_low, frequencies, low_floor, low_max)

    ax_high.scatter(
        modes,
        frequencies,
        c=bp,
        cmap=layer_cmap,
        norm=norm,
        **scatter_options,
    )
    ax_high.set_title(
        f"(d) Optical-mode zoom ({high_min:g}-{high_max:g} cm$^{{-1}}$)"
    )
    style_axis(ax_high, len(frequencies))
    ax_high.set_ylim(high_min, high_max)
    zoom_mode_axis_to_frequency_window(ax_high, frequencies, high_min, high_max)

    for ax in axes.ravel():
        mark_selected_modes(ax, modes, frequencies, selected_modes)

    fig.suptitle(
        r"$\Gamma$-point projected phonon modes of BP/MoS$_2$ (94 atoms)",
        fontsize=15,
        fontweight="semibold",
    )
    fig.text(
        0.5,
        -0.012,
        "Projection weights are normalized Phonopy eigenvector weights. "
        "This figure contains Gamma-point modes only and is not a q-dependent dispersion.",
        ha="center",
        va="top",
        fontsize=8.3,
        color="#444444",
    )

    png_path = output_prefix.with_suffix(".png")
    pdf_path = output_prefix.with_suffix(".pdf")
    fig.savefig(png_path, dpi=dpi, bbox_inches="tight", facecolor="white")
    fig.savefig(pdf_path, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    return png_path, pdf_path


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Plot BP/MoS2 layer and polarization projections for all "
            "Gamma-point modes in a Phonopy qpoints.yaml file."
        )
    )
    parser.add_argument(
        "qpoints",
        nargs="?",
        type=Path,
        default=Path("qpoints.yaml"),
        help="Phonopy qpoints.yaml file (default: qpoints.yaml)",
    )
    parser.add_argument(
        "--structure-yaml",
        type=Path,
        default=None,
        help=(
            "Phonopy structure YAML containing atom symbols; default is "
            "phonopy.yaml beside qpoints.yaml"
        ),
    )
    parser.add_argument(
        "--input-frequency-unit",
        choices=("thz", "cm-1"),
        default="thz",
        help=(
            "unit used by frequency values in qpoints.yaml "
            "(default: thz, Phonopy's standard output)"
        ),
    )
    parser.add_argument(
        "--output-prefix",
        type=Path,
        default=Path("gamma_phonon_projection"),
        help="output filename prefix (default: gamma_phonon_projection)",
    )
    parser.add_argument(
        "--label-modes",
        type=parse_mode_list,
        default=parse_mode_list("8,225"),
        help="comma-separated one-based modes to label (default: 8,225)",
    )
    parser.add_argument("--low-max", type=float, default=100.0)
    parser.add_argument("--high-min", type=float, default=330.0)
    parser.add_argument("--high-max", type=float, default=440.0)
    parser.add_argument("--dpi", type=int, default=300)
    return parser


def main() -> int:
    args = build_parser().parse_args()
    if args.high_min >= args.high_max:
        raise ValueError("--high-min must be smaller than --high-max")
    if not args.qpoints.is_file():
        raise FileNotFoundError(args.qpoints)

    args.output_prefix.parent.mkdir(parents=True, exist_ok=True)

    structure_yaml = (
        args.structure_yaml
        if args.structure_yaml is not None
        else args.qpoints.with_name("phonopy.yaml")
    )
    symbols, raw_frequencies, eigenvectors, gamma_q, symbol_source = load_gamma_modes(
        args.qpoints, structure_yaml
    )
    if args.input_frequency_unit == "thz":
        frequencies_thz = raw_frequencies
        frequencies = raw_frequencies * THZ_TO_CM1
    else:
        frequencies = raw_frequencies
        frequencies_thz = raw_frequencies / THZ_TO_CM1
    bp, mos2, x, y, z, in_plane = calculate_projections(symbols, eigenvectors)

    csv_path = args.output_prefix.with_suffix(".csv")
    write_csv(
        csv_path,
        frequencies,
        frequencies_thz,
        bp,
        mos2,
        x,
        y,
        z,
        in_plane,
    )
    png_path, pdf_path = plot_projection(
        args.output_prefix,
        frequencies,
        bp,
        z,
        args.label_modes,
        args.low_max,
        args.high_min,
        args.high_max,
        args.dpi,
    )

    counts = Counter(symbols.tolist())
    negative_modes = np.flatnonzero(frequencies < -0.01) + 1
    near_zero_modes = np.flatnonzero(np.abs(frequencies) <= 0.01) + 1
    print(
        f"Loaded Gamma q={gamma_q.tolist()}: {len(frequencies)} modes, "
        f"atoms={dict(counts)}; symbols from {symbol_source}"
    )
    if len(negative_modes):
        summary = ", ".join(
            f"{mode:03d} ({frequencies[mode - 1]:.4f} cm^-1)"
            for mode in negative_modes
        )
        print(f"Negative-frequency modes: {summary}")
    else:
        print("Negative-frequency modes: none")
    if len(near_zero_modes):
        summary = ", ".join(
            f"{mode:03d} ({frequencies[mode - 1]:.6f} cm^-1)"
            for mode in near_zero_modes
        )
        print(f"Numerically near-zero modes (|frequency| <= 0.01 cm^-1): {summary}")

    for mode in args.label_modes:
        if mode <= len(frequencies):
            index = mode - 1
            print(
                f"Mode {mode:03d}: {frequencies[index]:.4f} cm^-1 "
                f"({frequencies_thz[index]:.6f} THz), "
                f"BP={100 * bp[index]:.2f}%, MoS2={100 * mos2[index]:.2f}%, "
                f"z={100 * z[index]:.2f}%"
            )
        else:
            print(
                f"Warning: requested mode {mode} exceeds {len(frequencies)} modes",
                file=sys.stderr,
            )

    print(f"Wrote {png_path}")
    print(f"Wrote {pdf_path}")
    print(f"Wrote {csv_path}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (KeyError, TypeError, ValueError, OSError, yaml.YAMLError) as exc:
        print(f"Error: {exc}", file=sys.stderr)
        raise SystemExit(1)

