#!/usr/bin/env python3
"""
Compute mass-weighted percentiles for a PID-selected tracer population,
optionally applying cuts from stat_traj.dat.

The output files are directly readable by

    python3 plot_timescales.py --prefix PREFIX

Outputs
-------
  PREFIX_ye.txt
  PREFIX_yecap.txt
  PREFIX_texp.txt
  PREFIX_caprate.txt
  PREFIX_eta.txt
  PREFIX_r.txt
  PREFIX_vr.txt
  PREFIX_time.txt
  PREFIX_avtexp.txt
  PREFIX_ids.txt

The PID source may be either a one-column ID list or a wider *_particles.txt
table. In either case, the first non-comment column is interpreted as PID.

Final-state cuts
----------------
  --ye-min / --ye-max
  --s-min  / --s-max
  --vr-min / --vr-max

Alternative 5-GK cuts
-------------------
  --ye-5gk-min   / --ye-5gk-max
  --s-5gk-min    / --s-5gk-max
  --texp-5gk-min / --texp-5gk-max

Final-state cuts and 5-GK cuts are mutually exclusive.
"""

from __future__ import annotations

import argparse
from pathlib import Path
import sys

import numpy as np


DEFAULT_FILES = {
    "time": "traj_time",
    "texp": "traj_texp",
    "caprate": "traj_caprate",
    "ye": "traj_Ye",
    "yecap": "traj_yecap",
    "eta": "traj_eta",
    "r": "traj_r",
    "vr": "traj_vel",
    "avtexp": "traj_avtexp",
}


# stat_traj.dat columns, 0-based.
STAT_COL_ID = 0
STAT_COL_MASS = 1
STAT_COL_VR_FIN = 6
STAT_COL_S_FIN = 10
STAT_COL_YE_FIN = 11
STAT_COL_S_5GK = 13
STAT_COL_YE_5GK = 14
STAT_COL_TEXP_5GK = 15


def normalized_prefix(prefix: str | Path) -> Path:
    path = Path(prefix)
    if path.suffix.lower() in (".txt", ".dat"):
        path = path.with_suffix("")
    path.parent.mkdir(parents=True, exist_ok=True)
    return path


def output_filename(prefix: str | Path, quantity: str) -> Path:
    path = normalized_prefix(prefix)
    return path.parent / f"{path.name}_{quantity}.txt"


def ids_output_filename(prefix: str | Path) -> Path:
    path = normalized_prefix(prefix)
    return path.parent / f"{path.name}_ids.txt"


def plot_output_filename(
    prefix: str | Path,
    plot_name: str,
    extension: str,
) -> Path:
    path = normalized_prefix(prefix)
    return path.parent / f"{plot_name}-{path.name}.{extension.lstrip('.')}"


def read_pid_list(filename: str | Path) -> np.ndarray:
    """
    Read PIDs from the first non-comment column.

    Duplicate PIDs are removed while preserving their first occurrence.
    """
    filename = Path(filename)
    raw = np.loadtxt(
        filename,
        comments="#",
        usecols=0,
        ndmin=1,
    )

    if raw.size == 0:
        raise ValueError(f"{filename}: no particle IDs were found.")

    rounded = np.rint(raw)
    if not np.allclose(raw, rounded, rtol=0.0, atol=1.0e-8):
        raise ValueError(
            f"{filename}: the first column contains non-integer PIDs."
        )

    seen: set[int] = set()
    ids_unique: list[int] = []

    for pid in rounded.astype(np.int64):
        pid_int = int(pid)
        if pid_int not in seen:
            seen.add(pid_int)
            ids_unique.append(pid_int)

    duplicate_count = raw.size - len(ids_unique)
    if duplicate_count:
        print(
            f"Warning: removed {duplicate_count} duplicate PID entries "
            f"from {filename}.",
            file=sys.stderr,
        )

    return np.asarray(ids_unique, dtype=np.int64)


def read_temperature_table(
    filename: str | Path,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """
    Read a temperature-dependent tracer table.

    Format
    ------
      # optional comments
      T1 T2 T3 ...
      pid mass value(T1) value(T2) value(T3) ...
      ...
    """
    filename = Path(filename)
    rows: list[list[str]] = []

    with filename.open("r") as stream:
        for raw_line in stream:
            line = raw_line.strip()
            if not line or line.startswith("#"):
                continue
            rows.append(line.split())

    if len(rows) < 2:
        raise ValueError(
            f"{filename}: expected one temperature row and at least "
            "one tracer row."
        )

    temperatures = np.asarray(rows[0], dtype=float)
    ntemperature = temperatures.size

    if ntemperature == 0:
        raise ValueError(f"{filename}: empty temperature row.")

    ids: list[int] = []
    masses: list[float] = []
    values: list[list[float]] = []

    for line_number, row in enumerate(rows[1:], start=2):
        expected_columns = ntemperature + 2
        if len(row) != expected_columns:
            raise ValueError(
                f"{filename}: data line {line_number}: expected "
                f"{expected_columns} columns, got {len(row)}."
            )

        pid_float = float(row[0])
        pid_round = int(round(pid_float))
        if not np.isclose(pid_float, pid_round, rtol=0.0, atol=1.0e-8):
            raise ValueError(
                f"{filename}: data line {line_number}: PID is not integer."
            )

        ids.append(pid_round)
        masses.append(float(row[1]))
        values.append([float(item) for item in row[2:]])

    ids_array = np.asarray(ids, dtype=np.int64)
    if np.unique(ids_array).size != ids_array.size:
        raise ValueError(f"{filename}: duplicate tracer IDs are present.")

    return (
        temperatures,
        ids_array,
        np.asarray(masses, dtype=float),
        np.asarray(values, dtype=float),
    )


def read_stat_traj(
    filename: str | Path,
) -> dict[str, np.ndarray]:
    """Read the stat_traj.dat quantities used for cuts and weights."""
    filename = Path(filename)
    data = np.loadtxt(filename, comments="#", ndmin=2)

    ncol_needed = max(
        STAT_COL_ID,
        STAT_COL_MASS,
        STAT_COL_VR_FIN,
        STAT_COL_S_FIN,
        STAT_COL_YE_FIN,
        STAT_COL_S_5GK,
        STAT_COL_YE_5GK,
        STAT_COL_TEXP_5GK,
    ) + 1

    if data.shape[1] < ncol_needed:
        raise ValueError(
            f"{filename}: expected at least {ncol_needed} columns, "
            f"got {data.shape[1]}."
        )

    result = {
        "ids": data[:, STAT_COL_ID].astype(np.int64),
        "mass": data[:, STAT_COL_MASS].astype(float),
        "vr_fin": data[:, STAT_COL_VR_FIN].astype(float),
        "s_fin": data[:, STAT_COL_S_FIN].astype(float),
        "ye_fin": data[:, STAT_COL_YE_FIN].astype(float),
        "s_5gk": data[:, STAT_COL_S_5GK].astype(float),
        "ye_5gk": data[:, STAT_COL_YE_5GK].astype(float),
        "texp_5gk": data[:, STAT_COL_TEXP_5GK].astype(float),
    }

    if np.unique(result["ids"]).size != result["ids"].size:
        raise ValueError(f"{filename}: duplicate tracer IDs are present.")

    return result


def make_pid_index(
    ids: np.ndarray,
    *,
    label: str,
) -> dict[int, int]:
    ids = np.asarray(ids, dtype=np.int64)

    if np.unique(ids).size != ids.size:
        raise ValueError(f"{label}: duplicate tracer IDs are present.")

    return {int(pid): index for index, pid in enumerate(ids)}


def report_missing_ids(
    requested_ids: np.ndarray,
    available_ids: np.ndarray,
    *,
    label: str,
    mode: str,
) -> list[int]:
    available = set(int(pid) for pid in available_ids)
    missing = [
        int(pid)
        for pid in requested_ids
        if int(pid) not in available
    ]

    if not missing:
        return []

    message = (
        f"{label}: {len(missing)} requested PIDs are missing. "
        f"First missing PID: {missing[0]}"
    )

    if mode == "error":
        raise ValueError(message)
    if mode == "warn":
        print(f"Warning: {message}", file=sys.stderr)

    return missing


def common_selected_ids(
    requested_ids: np.ndarray,
    available_id_sets: list[tuple[str, np.ndarray]],
    *,
    missing_mode: str,
) -> np.ndarray:
    """Use the PID intersection of all required inputs."""
    common = set(int(pid) for pid in requested_ids)

    for label, available_ids in available_id_sets:
        report_missing_ids(
            requested_ids,
            available_ids,
            label=label,
            mode=missing_mode,
        )
        common &= set(int(pid) for pid in available_ids)

    selected = np.asarray(
        [int(pid) for pid in requested_ids if int(pid) in common],
        dtype=np.int64,
    )

    if selected.size == 0:
        raise ValueError(
            "No selected PID is present in every required input file."
        )

    return selected


def reorder_rows(
    requested_ids: np.ndarray,
    source_ids: np.ndarray,
    *source_arrays: np.ndarray,
    label: str,
) -> tuple[np.ndarray, ...]:
    """Reorder source arrays into requested PID order."""
    source_index = make_pid_index(source_ids, label=label)

    missing = [
        int(pid)
        for pid in requested_ids
        if int(pid) not in source_index
    ]
    if missing:
        raise ValueError(
            f"{label}: {len(missing)} selected PIDs are missing. "
            f"First missing PID: {missing[0]}"
        )

    order = np.fromiter(
        (source_index[int(pid)] for pid in requested_ids),
        dtype=np.int64,
        count=requested_ids.size,
    )

    return tuple(array[order] for array in source_arrays)


def add_range_condition(
    mask: np.ndarray,
    values: np.ndarray,
    minimum: float | None,
    maximum: float | None,
    *,
    open_upper: bool,
) -> np.ndarray:
    """Apply a finite-value and optional lower/upper cut."""
    condition = np.isfinite(values)

    if minimum is not None:
        condition &= values >= minimum

    if maximum is not None:
        if open_upper:
            condition &= values < maximum
        else:
            condition &= values <= maximum

    return mask & condition


def validate_range(
    name: str,
    minimum: float | None,
    maximum: float | None,
    *,
    logarithmic: bool = False,
) -> None:
    for suffix, value in (("min", minimum), ("max", maximum)):
        if value is not None and not np.isfinite(value):
            raise ValueError(f"{name}_{suffix} must be finite.")

        if logarithmic and value is not None and value <= 0.0:
            raise ValueError(
                f"{name}_{suffix} must be positive for a logarithmic quantity."
            )

    if (
        minimum is not None
        and maximum is not None
        and minimum >= maximum
    ):
        raise ValueError(
            f"{name}_min must be smaller than {name}_max."
        )


def weighted_quantile(
    values: np.ndarray,
    weights: np.ndarray,
    quantiles: np.ndarray,
) -> np.ndarray:
    """Mass-weighted quantiles using the weighted CDF at tracer-bin centers."""
    values = np.asarray(values, dtype=float)
    weights = np.asarray(weights, dtype=float)
    quantiles = np.asarray(quantiles, dtype=float)

    valid = (
        np.isfinite(values)
        & np.isfinite(weights)
        & (weights > 0.0)
    )
    values = values[valid]
    weights = weights[valid]

    if values.size == 0 or np.sum(weights) <= 0.0:
        return np.full(quantiles.shape, np.nan, dtype=float)

    order = np.argsort(values)
    values_sorted = values[order]
    weights_sorted = weights[order]

    cdf = (
        np.cumsum(weights_sorted) - 0.5 * weights_sorted
    ) / np.sum(weights_sorted)

    return np.interp(
        quantiles,
        cdf,
        values_sorted,
        left=values_sorted[0],
        right=values_sorted[-1],
    )


def percentile_label(percentile: float) -> str:
    return "p" + f"{percentile:g}".replace(".", "p")


def make_aligned_header(
    header_columns: list[str],
    widths: list[int],
    comment: str = "# ",
) -> str:
    parts = []

    for index, (name, width) in enumerate(
        zip(header_columns, widths)
    ):
        if index == 0:
            parts.append(f"{name:>{width - len(comment)}s}")
        else:
            parts.append(f"{name:>{width}s}")

    return comment + " ".join(parts)


def write_percentile_table(
    output_file: str | Path,
    *,
    quantity_name: str,
    temperatures: np.ndarray,
    values: np.ndarray,
    weights: np.ndarray,
    percentiles: np.ndarray,
) -> None:
    """
    Write the column layout used by calc_percentiles_by_stat_cuts.py.

    With the default percentiles, plot_timescales.py reads:
      col 0: temperature
      col 1: lower percentile
      col 2: median
      col 3: upper percentile
    """
    quantiles = percentiles / 100.0
    rows = []

    for column, temperature in enumerate(temperatures):
        values_at_temperature = values[:, column]
        valid = (
            np.isfinite(values_at_temperature)
            & np.isfinite(weights)
            & (weights > 0.0)
        )

        quantile_values = weighted_quantile(
            values_at_temperature[valid],
            weights[valid],
            quantiles,
        )

        selected_mass = np.sum(weights[valid])
        n_selected = np.count_nonzero(valid)

        rows.append(
            np.concatenate(
                (
                    [temperature],
                    quantile_values,
                    [selected_mass, n_selected],
                )
            )
        )

    rows_array = np.asarray(rows, dtype=float)

    quantile_labels = [
        f"{quantity_name}_{percentile_label(p)}_mass_weighted"
        for p in percentiles
    ]
    header_columns = (
        ["T"]
        + quantile_labels
        + ["selected_mass", "n_selected"]
    )

    widths = [22]
    widths += [
        max(28, len(label) + 2)
        for label in quantile_labels
    ]
    widths += [22, 12]

    header = make_aligned_header(
        header_columns,
        widths,
        comment="# ",
    )

    formats = " ".join(
        [f"%{widths[0]}.8e"]
        + [
            f"%{widths[index]}.8e"
            for index in range(1, 1 + percentiles.size)
        ]
        + [
            f"%{widths[-2]}.8e",
            f"%{widths[-1]}d",
        ]
    )

    np.savetxt(
        output_file,
        rows_array,
        header=header,
        fmt=formats,
        comments="",
    )
    print(f"Wrote: {output_file}")


def write_id_file(
    output_file: str | Path,
    ids: np.ndarray,
) -> None:
    np.savetxt(
        output_file,
        np.asarray(ids, dtype=np.int64),
        fmt="%d",
        header="IDs for selected trajectories",
    )
    print(f"Wrote: {output_file}")


def make_selection_plots(
    *,
    prefix: str | Path,
    plot_extension: str,
    cut_mode: str,
    all_ye_fin: np.ndarray,
    all_s_fin: np.ndarray,
    all_vr_fin: np.ndarray,
    all_ye_5gk: np.ndarray,
    all_s_5gk: np.ndarray,
    all_texp_5gk: np.ndarray,
    pid_ye_fin: np.ndarray,
    pid_s_fin: np.ndarray,
    pid_vr_fin: np.ndarray,
    pid_ye_5gk: np.ndarray,
    pid_s_5gk: np.ndarray,
    pid_texp_5gk: np.ndarray,
    ye_min: float | None,
    ye_max: float | None,
    s_min: float | None,
    s_max: float | None,
    vr_min: float | None,
    vr_max: float | None,
    ye_5gk_min: float | None,
    ye_5gk_max: float | None,
    s_5gk_min: float | None,
    s_5gk_max: float | None,
    texp_5gk_min: float | None,
    texp_5gk_max: float | None,
) -> None:
    """Draw the active cut-space scatter plots."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.ticker as ticker
    from matplotlib.patches import Rectangle

    params={
        'font.size'           : 24.0     ,
        #    'font.family'         : 'DeJaVu Sans'  ,
        'font.family'         : 'Times New Roman'  ,
        'axes.grid.axis'      : 'both',
        'xtick.major.size'    : 6        ,
        'xtick.major.width'   : 1.5      ,
        'xtick.minor.size'    : 3        ,
        'xtick.minor.width'   : 1.5      ,
        'xtick.labelsize'     : 24.0     ,
        'xtick.direction'     : 'out'     ,
        'ytick.major.size'    : 6        ,
        'ytick.major.width'   : 1.5      ,
        'ytick.minor.size'    : 3        ,
        'ytick.minor.width'   : 1.5      ,
        'ytick.labelsize'     : 24.0     ,
        'ytick.direction'     : 'out'     ,
        'xtick.major.pad'     : 2        ,
        'xtick.minor.pad'     : 2        ,
        'ytick.major.pad'     : 2        ,
        'ytick.minor.pad'     : 2        ,
        'axes.linewidth'      : 1.5      ,
        'text.usetex'         : True    }
    
    plt.rcParams.update(params)


    if cut_mode not in ("final", "5gk"):
        raise ValueError(f"Unknown cut mode: {cut_mode}")

    def draw_plane(
        all_x: np.ndarray,
        all_y: np.ndarray,
        pid_x: np.ndarray,
        pid_y: np.ndarray,
        *,
        xlabel: str,
        ylabel: str,
        x_display: tuple[float, float],
        y_display: tuple[float, float],
        x_cut: tuple[float | None, float | None],
        y_cut: tuple[float | None, float | None],
        xscale: str,
        yscale: str,
        output_name: str,
    ) -> None:
        all_finite = np.isfinite(all_x) & np.isfinite(all_y)
        pid_finite = np.isfinite(pid_x) & np.isfinite(pid_y)

        if xscale == "log":
            all_finite &= all_x > 0.0
            pid_finite &= pid_x > 0.0
        if yscale == "log":
            all_finite &= all_y > 0.0
            pid_finite &= pid_y > 0.0

        fig = plt.figure(figsize=(10.0, 6.25))
        ax = fig.add_subplot(111)

        ax.scatter(
            all_x[all_finite],
            all_y[all_finite],
            marker=".",
            s=3,
            color="0.75",
            linewidths=0.0,
            alpha=0.8,
            label="all stat_traj.dat",
            zorder=1,
        )

        ax.scatter(
            pid_x[pid_finite],
            pid_y[pid_finite],
            marker=".",
            s=5,
            facecolors="none",
            edgecolors="black",
            linewidths=0.7,
            alpha=0.9,
            label="input PID list",
            zorder=3,
        )

        x0 = x_cut[0] if x_cut[0] is not None else x_display[0]
        x1 = x_cut[1] if x_cut[1] is not None else x_display[1]
        y0 = y_cut[0] if y_cut[0] is not None else y_display[0]
        y1 = y_cut[1] if y_cut[1] is not None else y_display[1]

        if x0 < x1 and y0 < y1:
            rectangle = Rectangle(
                (x0, y0),
                x1 - x0,
                y1 - y0,
                fill=False,
                edgecolor="red",
                linewidth=2.0,
                label="cut range",
                zorder=4,
            )
            ax.add_patch(rectangle)

        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.set_xscale(xscale)
        ax.set_yscale(yscale)
        ax.set_xlim(*x_display)
        ax.set_ylim(*y_display)

        if xscale == "log":
            ax.xaxis.set_major_locator(ticker.LogLocator(numticks=15))
            ax.xaxis.set_minor_locator(
                ticker.LogLocator(
                    numticks=15,
                    subs=(0.1, 0.2, 0.3, 0.4, 0.5,
                          0.6, 0.7, 0.8, 0.9),
                )
            )
        if yscale == "log":
            ax.yaxis.set_major_locator(ticker.LogLocator(numticks=15))
            ax.yaxis.set_minor_locator(
                ticker.LogLocator(
                    numticks=15,
                    subs=(0.1, 0.2, 0.3, 0.4, 0.5,
                          0.6, 0.7, 0.8, 0.9),
                )
            )

        ax.tick_params(axis="both", which="major", direction="out", width=1.5, length=6)
        ax.tick_params(axis="both", which="minor", direction="out", width=1.5, length=3)
        ax.legend(fontsize=11, frameon=False, loc="best")

        fig.subplots_adjust(left=0.15, bottom=0.15, right=0.95, top=0.95)
        output_file = plot_output_filename(prefix, output_name, plot_extension)
        fig.savefig(output_file)
        plt.close(fig)
        print(f"Wrote: {output_file}")

    if cut_mode == "final":
        draw_plane(
            all_ye_fin, all_s_fin, pid_ye_fin, pid_s_fin,
            xlabel=r"$Y_{\rm e,fin}$",
            ylabel=r"$s_{\rm fin}/k_{\rm B}$",
            x_display=(0.0, 0.6), y_display=(1.0, 3.0e2),
            x_cut=(ye_min, ye_max), y_cut=(s_min, s_max),
            xscale="linear", yscale="log", output_name="Ye-S",
        )
        draw_plane(
            all_ye_fin, all_vr_fin, pid_ye_fin, pid_vr_fin,
            xlabel=r"$Y_{\rm e,fin}$",
            ylabel=r"$v^r_{\rm fin}\ [{\rm cm\,s^{-1}}]$",
            x_display=(0.0, 0.6), y_display=(1.0e8, 3.0e10),
            x_cut=(ye_min, ye_max), y_cut=(vr_min, vr_max),
            xscale="linear", yscale="log", output_name="Ye-vr",
        )
        draw_plane(
            all_s_fin, all_vr_fin, pid_s_fin, pid_vr_fin,
            xlabel=r"$s_{\rm fin}/k_{\rm B}$",
            ylabel=r"$v^r_{\rm fin}\ [{\rm cm\,s^{-1}}]$",
            x_display=(1.0, 3.0e2), y_display=(1.0e8, 3.0e10),
            x_cut=(s_min, s_max), y_cut=(vr_min, vr_max),
            xscale="log", yscale="log", output_name="S-vr",
        )
    else:
        draw_plane(
            all_ye_5gk, all_s_5gk, pid_ye_5gk, pid_s_5gk,
            xlabel=r"$Y_{\rm e}(5\,{\rm GK})$",
            ylabel=r"$s(5\,{\rm GK})/k_{\rm B}$",
            x_display=(0.0, 0.6), y_display=(1.0, 3.0e2),
            x_cut=(ye_5gk_min, ye_5gk_max), y_cut=(s_5gk_min, s_5gk_max),
            xscale="linear", yscale="log", output_name="Ye-S_5GK",
        )
        draw_plane(
            all_ye_5gk, all_texp_5gk, pid_ye_5gk, pid_texp_5gk,
            xlabel=r"$Y_{\rm e}(5\,{\rm GK})$",
            ylabel=r"$t_{\rm exp}(5\,{\rm GK})\ [{\rm s}]$",
            x_display=(0.0, 0.6), y_display=(1.0e-4, 1.0e1),
            x_cut=(ye_5gk_min, ye_5gk_max), y_cut=(texp_5gk_min, texp_5gk_max),
            xscale="linear", yscale="log", output_name="Ye-texp_5GK",
        )
        draw_plane(
            all_s_5gk, all_texp_5gk, pid_s_5gk, pid_texp_5gk,
            xlabel=r"$s(5\,{\rm GK})/k_{\rm B}$",
            ylabel=r"$t_{\rm exp}(5\,{\rm GK})\ [{\rm s}]$",
            x_display=(1.0, 3.0e2), y_display=(1.0e-4, 1.0e1),
            x_cut=(s_5gk_min, s_5gk_max), y_cut=(texp_5gk_min, texp_5gk_max),
            xscale="log", yscale="log", output_name="S-texp_5GK",
        )


def add_cut_arguments(parser: argparse.ArgumentParser) -> None:
    # Final-state cuts.
    parser.add_argument(
        "--ye-min",
        "--ye_min",
        dest="ye_min",
        type=float,
        default=None,
        help="Minimum Ye_fin. Omit for no lower cut.",
    )
    parser.add_argument(
        "--ye-max",
        "--ye_max",
        dest="ye_max",
        type=float,
        default=None,
        help="Maximum Ye_fin. Omit for no upper cut.",
    )
    parser.add_argument(
        "--s-min",
        "--s_min",
        dest="s_min",
        type=float,
        default=None,
        help="Minimum s_fin [k_B/nuc]. Omit for no lower cut.",
    )
    parser.add_argument(
        "--s-max",
        "--s_max",
        dest="s_max",
        type=float,
        default=None,
        help="Maximum s_fin [k_B/nuc]. Omit for no upper cut.",
    )
    parser.add_argument(
        "--vr-min",
        "--vr_min",
        dest="vr_min",
        type=float,
        default=None,
        help="Minimum v^r_fin [cm/s]. Omit for no lower cut.",
    )
    parser.add_argument(
        "--vr-max",
        "--vr_max",
        dest="vr_max",
        type=float,
        default=None,
        help="Maximum v^r_fin [cm/s]. Omit for no upper cut.",
    )

    # Alternative 5-GK cuts; do not combine with final-state cuts.
    parser.add_argument(
        "--ye-5gk-min",
        "--ye_5gk_min",
        "--ye-5GK-min",
        dest="ye_5gk_min",
        type=float,
        default=None,
        help="Minimum Ye at 5 GK. Omit for no lower cut.",
    )
    parser.add_argument(
        "--ye-5gk-max",
        "--ye_5gk_max",
        "--ye-5GK-max",
        dest="ye_5gk_max",
        type=float,
        default=None,
        help="Maximum Ye at 5 GK. Omit for no upper cut.",
    )
    parser.add_argument(
        "--s-5gk-min",
        "--s_5gk_min",
        "--s-5GK-min",
        dest="s_5gk_min",
        type=float,
        default=None,
        help="Minimum entropy at 5 GK [k_B/nuc].",
    )
    parser.add_argument(
        "--s-5gk-max",
        "--s_5gk_max",
        "--s-5GK-max",
        dest="s_5gk_max",
        type=float,
        default=None,
        help="Maximum entropy at 5 GK [k_B/nuc].",
    )
    parser.add_argument(
        "--texp-5gk-min",
        "--texp_5gk_min",
        "--t-exp-5gk-min",
        "--t_exp_5gk_min",
        dest="texp_5gk_min",
        type=float,
        default=None,
        help="Minimum expansion time at 5 GK [s].",
    )
    parser.add_argument(
        "--texp-5gk-max",
        "--texp_5gk_max",
        "--t-exp-5gk-max",
        "--t_exp_5gk_max",
        dest="texp_5gk_max",
        type=float,
        default=None,
        help="Maximum expansion time at 5 GK [s].",
    )

    parser.add_argument(
        "--open-upper",
        "--open_upper",
        dest="open_upper",
        action="store_true",
        help="Use < rather than <= for every upper bound.",
    )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Compute mass-weighted percentile tables for a PID-selected "
            "population after optional final-state and 5-GK cuts."
        )
    )

    parser.add_argument(
        "--id-file",
        "--id_file",
        "--pid-file",
        "--pid_file",
        "--particle-file",
        "--particle_file",
        dest="id_file",
        required=True,
        help=(
            "PID source file. The first non-comment column is used, "
            "so a *_particles.txt table can be supplied directly."
        ),
    )
    parser.add_argument(
        "--stat-file",
        "--stat_file",
        dest="stat_file",
        default="stat_traj.dat",
        help="stat_traj.dat. Default: stat_traj.dat",
    )

    for quantity_name, default_filename in DEFAULT_FILES.items():
        parser.add_argument(
            f"--{quantity_name}-file",
            f"--{quantity_name}_file",
            dest=f"{quantity_name}_file",
            default=default_filename,
            help=(
                f"Input table for {quantity_name}. "
                f"Default: {default_filename}"
            ),
        )

    add_cut_arguments(parser)

    parser.add_argument(
        "--quantities",
        nargs="+",
        choices=tuple(DEFAULT_FILES),
        default=list(DEFAULT_FILES),
        help=(
            "Quantities to process. Default: all quantities needed by "
            "plot_timescales.py."
        ),
    )
    parser.add_argument(
        "--percentiles",
        type=float,
        nargs="+",
        default=[15.865, 50.0, 84.135],
        help=(
            "Mass-weighted percentiles in percent. "
            "Default: 15.865 50 84.135."
        ),
    )
    parser.add_argument(
        "--weight-source",
        "--weight_source",
        dest="weight_source",
        choices=("stat", "value"),
        default="stat",
        help=(
            "Use masses from stat_traj.dat or from each value table. "
            "Default: stat."
        ),
    )
    parser.add_argument(
        "--missing-ids",
        "--missing_ids",
        dest="missing_ids",
        choices=("error", "warn", "ignore"),
        default="error",
        help=(
            "Behavior when a requested PID is absent from an input. "
            "Default: error."
        ),
    )
    parser.add_argument(
        "--mass-rtol",
        "--mass_rtol",
        dest="mass_rtol",
        type=float,
        default=1.0e-10,
        help="Relative tolerance for comparing stat and table masses.",
    )
    parser.add_argument(
        "--plot-format",
        "--plot_format",
        dest="plot_format",
        default="png",
        help="Scatter-plot extension. Default: png.",
    )
    parser.add_argument(
        "--no-plots",
        "--no_plots",
        dest="no_plots",
        action="store_true",
        help="Do not write the selection scatter plots.",
    )
    parser.add_argument(
        "-o",
        "--output-prefix",
        "--output_prefix",
        dest="output_prefix",
        default="percentiles",
        help=(
            "Output prefix. Writes <prefix>_<quantity>.txt, "
            "<prefix>_ids.txt, and selection plots."
        ),
    )

    return parser.parse_args()


def main() -> None:
    args = parse_args()

    percentiles = np.asarray(args.percentiles, dtype=float)

    if percentiles.ndim != 1 or percentiles.size == 0:
        raise ValueError(
            "--percentiles must contain at least one value."
        )
    if np.any(~np.isfinite(percentiles)):
        raise ValueError("All percentiles must be finite.")
    if np.any((percentiles < 0.0) | (percentiles > 100.0)):
        raise ValueError(
            "All percentiles must be between 0 and 100."
        )

    if percentiles.size != 3:
        print(
            "Warning: plot_timescales.py assumes exactly three percentile "
            "columns after temperature. Direct use is safest with the "
            "default 15.865 50 84.135 values.",
            file=sys.stderr,
        )

    validate_range("ye", args.ye_min, args.ye_max)
    validate_range("s", args.s_min, args.s_max, logarithmic=True)
    validate_range("vr", args.vr_min, args.vr_max, logarithmic=True)
    validate_range(
        "ye_5gk",
        args.ye_5gk_min,
        args.ye_5gk_max,
    )
    validate_range(
        "s_5gk",
        args.s_5gk_min,
        args.s_5gk_max,
        logarithmic=True,
    )
    validate_range(
        "texp_5gk",
        args.texp_5gk_min,
        args.texp_5gk_max,
        logarithmic=True,
    )

    final_cut_requested = any(
        value is not None
        for value in (
            args.ye_min, args.ye_max,
            args.s_min, args.s_max,
            args.vr_min, args.vr_max,
        )
    )
    gk5_cut_requested = any(
        value is not None
        for value in (
            args.ye_5gk_min, args.ye_5gk_max,
            args.s_5gk_min, args.s_5gk_max,
            args.texp_5gk_min, args.texp_5gk_max,
        )
    )

    if final_cut_requested and gk5_cut_requested:
        raise ValueError(
            "Final-state cuts and 5-GK-state cuts cannot be used "
            "simultaneously. Specify only one cut family."
        )

    cut_mode = "5gk" if gk5_cut_requested else "final"
    print(f"Cut mode: {cut_mode}")

    requested_ids = read_pid_list(args.id_file)
    print(f"PID-list population: n={requested_ids.size}")

    stat = read_stat_traj(args.stat_file)

    # First reduce the PID source to rows available in stat_traj.dat.
    stat_available_ids = common_selected_ids(
        requested_ids,
        [(args.stat_file, stat["ids"])],
        missing_mode=args.missing_ids,
    )

    (
        stat_mass,
        ye_fin,
        s_fin,
        vr_fin,
        ye_5gk,
        s_5gk,
        texp_5gk,
    ) = reorder_rows(
        stat_available_ids,
        stat["ids"],
        stat["mass"],
        stat["ye_fin"],
        stat["s_fin"],
        stat["vr_fin"],
        stat["ye_5gk"],
        stat["s_5gk"],
        stat["texp_5gk"],
        label=args.stat_file,
    )

    cut_mask = (
        np.isfinite(stat_mass)
        & (stat_mass > 0.0)
    )

    if cut_mode == "final":
        cut_mask = add_range_condition(
            cut_mask, ye_fin, args.ye_min, args.ye_max,
            open_upper=args.open_upper,
        )
        cut_mask = add_range_condition(
            cut_mask, s_fin, args.s_min, args.s_max,
            open_upper=args.open_upper,
        )
        cut_mask = add_range_condition(
            cut_mask, vr_fin, args.vr_min, args.vr_max,
            open_upper=args.open_upper,
        )
    else:
        cut_mask = add_range_condition(
            cut_mask, ye_5gk, args.ye_5gk_min, args.ye_5gk_max,
            open_upper=args.open_upper,
        )
        cut_mask = add_range_condition(
            cut_mask, s_5gk, args.s_5gk_min, args.s_5gk_max,
            open_upper=args.open_upper,
        )
        cut_mask = add_range_condition(
            cut_mask, texp_5gk, args.texp_5gk_min, args.texp_5gk_max,
            open_upper=args.open_upper,
        )

    cut_selected_ids = stat_available_ids[cut_mask]

    if cut_selected_ids.size == 0:
        raise ValueError("No tracer remains after the requested cuts.")

    print(
        "After stat cuts: "
        f"n={cut_selected_ids.size}, "
        f"mass={np.sum(stat_mass[cut_mask]):.10e} g"
    )

    # Read all requested temperature-dependent tables and enforce a common
    # temperature grid, as required by plot_timescales.py.
    tables: dict[
        str,
        tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray],
    ] = {}

    available_id_sets: list[tuple[str, np.ndarray]] = []
    reference_temperatures: np.ndarray | None = None
    reference_temperature_label: str | None = None

    for quantity_name in args.quantities:
        filename = getattr(args, f"{quantity_name}_file")
        table = read_temperature_table(filename)
        temperatures, table_ids, _, _ = table

        if reference_temperatures is None:
            reference_temperatures = temperatures
            reference_temperature_label = filename
        else:
            same_grid = (
                temperatures.shape == reference_temperatures.shape
                and np.allclose(
                    temperatures,
                    reference_temperatures,
                    rtol=1.0e-12,
                    atol=0.0,
                    equal_nan=True,
                )
            )
            if not same_grid:
                raise ValueError(
                    f"{filename}: temperature grid differs from "
                    f"{reference_temperature_label}. "
                    "plot_timescales.py requires a common grid."
                )

        tables[quantity_name] = table
        available_id_sets.append((filename, table_ids))

    final_ids = common_selected_ids(
        cut_selected_ids,
        available_id_sets,
        missing_mode=args.missing_ids,
    )

    print(f"Final common selected population: n={final_ids.size}")

    write_id_file(
        ids_output_filename(args.output_prefix),
        final_ids,
    )

    if not args.no_plots:
        make_selection_plots(
            prefix=args.output_prefix,
            plot_extension=args.plot_format,
            cut_mode=cut_mode,
            all_ye_fin=stat["ye_fin"],
            all_s_fin=stat["s_fin"],
            all_vr_fin=stat["vr_fin"],
            all_ye_5gk=stat["ye_5gk"],
            all_s_5gk=stat["s_5gk"],
            all_texp_5gk=stat["texp_5gk"],
            pid_ye_fin=ye_fin,
            pid_s_fin=s_fin,
            pid_vr_fin=vr_fin,
            pid_ye_5gk=ye_5gk,
            pid_s_5gk=s_5gk,
            pid_texp_5gk=texp_5gk,
            ye_min=args.ye_min,
            ye_max=args.ye_max,
            s_min=args.s_min,
            s_max=args.s_max,
            vr_min=args.vr_min,
            vr_max=args.vr_max,
            ye_5gk_min=args.ye_5gk_min,
            ye_5gk_max=args.ye_5gk_max,
            s_5gk_min=args.s_5gk_min,
            s_5gk_max=args.s_5gk_max,
            texp_5gk_min=args.texp_5gk_min,
            texp_5gk_max=args.texp_5gk_max,
        )

    # Stat weights in exact final-ID order.
    (final_stat_weights,) = reorder_rows(
        final_ids,
        stat["ids"],
        stat["mass"],
        label=args.stat_file,
    )

    for quantity_name in args.quantities:
        filename = getattr(args, f"{quantity_name}_file")
        (
            temperatures,
            table_ids,
            table_masses,
            table_values,
        ) = tables[quantity_name]

        (
            selected_table_masses,
            selected_values,
        ) = reorder_rows(
            final_ids,
            table_ids,
            table_masses,
            table_values,
            label=filename,
        )

        if args.weight_source == "stat":
            weights = final_stat_weights

            valid_comparison = (
                np.isfinite(weights)
                & np.isfinite(selected_table_masses)
            )
            if np.any(valid_comparison):
                denominator = np.maximum(
                    np.abs(weights[valid_comparison]),
                    1.0,
                )
                max_relative_difference = np.max(
                    np.abs(
                        selected_table_masses[valid_comparison]
                        - weights[valid_comparison]
                    )
                    / denominator
                )

                if max_relative_difference > args.mass_rtol:
                    print(
                        f"Warning: {filename}: table and stat masses differ; "
                        f"max relative difference="
                        f"{max_relative_difference:.3e}.",
                        file=sys.stderr,
                    )
        else:
            weights = selected_table_masses

        write_percentile_table(
            output_filename(
                args.output_prefix,
                quantity_name,
            ),
            quantity_name=quantity_name,
            temperatures=temperatures,
            values=selected_values,
            weights=weights,
            percentiles=percentiles,
        )

        print(
            f"  {quantity_name}: "
            f"n_selected={final_ids.size}, "
            f"n_temperature={temperatures.size}"
        )

    upper_symbol = "<" if args.open_upper else "<="
    print("Applied cuts:")

    if cut_mode == "final":
        print(
            f"  Ye_fin : {args.ye_min} <= Ye_fin "
            f"{upper_symbol} {args.ye_max}"
        )
        print(
            f"  s_fin  : {args.s_min} <= s_fin "
            f"{upper_symbol} {args.s_max}"
        )
        print(
            f"  vr_fin : {args.vr_min} <= vr_fin "
            f"{upper_symbol} {args.vr_max}"
        )
    else:
        print(
            f"  Ye_5GK   : {args.ye_5gk_min} <= Ye_5GK "
            f"{upper_symbol} {args.ye_5gk_max}"
        )
        print(
            f"  s_5GK    : {args.s_5gk_min} <= s_5GK "
            f"{upper_symbol} {args.s_5gk_max}"
        )
        print(
            f"  texp_5GK : {args.texp_5gk_min} <= texp_5GK "
            f"{upper_symbol} {args.texp_5gk_max}"
        )



if __name__ == "__main__":
    main()
