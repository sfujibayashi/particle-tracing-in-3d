#!/usr/bin/env python3
import argparse
from pathlib import Path
import numpy as np
import sys

import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as colors
import matplotlib.ticker as ticker
from matplotlib.patches import Rectangle

# stat_traj.dat columns, 0-based indices.
# Header in the uploaded file:
#  1 id
#  2 mass [g]
#  7 v^r_fin [cm/s]
# 11 s_fin [k_b/nuc]
# 12 Ye_fin
STAT_COL_ID = 0
STAT_COL_MASS = 1
STAT_COL_VR_FIN = 6
STAT_COL_S_FIN = 10
STAT_COL_YE_FIN = 11

STAT_COL_S_5GK = 13
STAT_COL_YE_5GK = 14
STAT_COL_HUT = 18
STAT_COL_X = 3
STAT_COL_Y = 4
STAT_COL_Z = 5
STAT_COL_T = 2

def read_temperature_table(filename):
    """
    Read a table with the following format:

      # id mass T
      T1 T2 T3 ...
      tracer_id tracer_mass value(T1) value(T2) value(T3) ...
      ...

    The values can include NaN.

    Returns
    -------
    temps : ndarray, shape (ntemp,)
    ids   : ndarray, shape (ntracer,)
    mass  : ndarray, shape (ntracer,)
    vals  : ndarray, shape (ntracer, ntemp)
    """
    rows = []
    with open(filename, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith("#"):
                continue
            rows.append(line.split())

    if len(rows) < 2:
        raise ValueError(f"{filename}: file must contain a temperature row and tracer rows.")

    temps = np.array(rows[0], dtype=float)
    ntemp = temps.size

    ids = []
    masses = []
    vals = []

    for iline, row in enumerate(rows[1:], start=2):
        if len(row) != ntemp + 2:
            raise ValueError(
                f"{filename}: line {iline}: expected {ntemp + 2} columns "
                f"(id, mass, {ntemp} values), got {len(row)}."
            )

        ids.append(int(float(row[0])))
        masses.append(float(row[1]))
        vals.append([float(x) for x in row[2:]])

    return (
        temps,
        np.array(ids, dtype=np.int64),
        np.array(masses, dtype=float),
        np.array(vals, dtype=float),
    )


def read_stat_traj(filename):
    """
    Read stat_traj.dat and extract id, mass, Ye_fin, s_fin, v^r_fin.

    Expected columns are based on the header:
      col  1: id
      col  2: mass [g]
      col  7: v^r_fin [cm/s]
      col 11: s_fin [k_b/nuc]
      col 12: Ye_fin
    """
    arr = np.loadtxt(filename, comments="#")

    if arr.ndim == 1:
        arr = arr[None, :]

    ncol_needed = max(STAT_COL_ID, STAT_COL_MASS, STAT_COL_VR_FIN, STAT_COL_S_FIN, STAT_COL_YE_FIN, STAT_COL_S_5GK, STAT_COL_YE_5GK, STAT_COL_HUT, STAT_COL_X, STAT_COL_Y, STAT_COL_Z, STAT_COL_T) + 1
    if arr.shape[1] < ncol_needed:
        raise ValueError(
            f"{filename}: expected at least {ncol_needed} columns, got {arr.shape[1]}."
        )

    ids = arr[:, STAT_COL_ID].astype(np.int64)
    mass = arr[:, STAT_COL_MASS].astype(float)
    vr_fin = arr[:, STAT_COL_VR_FIN].astype(float)
    s_fin = arr[:, STAT_COL_S_FIN].astype(float)
    ye_fin = arr[:, STAT_COL_YE_FIN].astype(float)
    s_5GK = arr[:, STAT_COL_S_5GK].astype(float)
    ye_5GK = arr[:, STAT_COL_YE_5GK].astype(float)
    bernoulli = arr[:, STAT_COL_HUT].astype(float)
    x_fin = arr[:, STAT_COL_X].astype(float)
    y_fin = arr[:, STAT_COL_Y].astype(float)
    z_fin = arr[:, STAT_COL_Z].astype(float)
    t_fin = arr[:, STAT_COL_T].astype(float)

    r_fin = np.sqrt(x_fin**2 + y_fin**2 + z_fin**2)

    return ids, mass, ye_fin, s_fin, vr_fin, ye_5GK, s_5GK, bernoulli, r_fin, t_fin


def align_to_reference_ids(ids_ref, ids_other, *arrays_other, label="other"):
    """
    Reorder arrays_other so that their first axis matches ids_ref.

    Returns reordered arrays in the same order as arrays_other.
    """
    other_index = {tid: i for i, tid in enumerate(ids_other)}

    missing = [tid for tid in ids_ref if tid not in other_index]
    if missing:
        raise ValueError(
            f"{len(missing)} tracer IDs in the reference file are missing from {label}. "
            f"First missing ID: {missing[0]}"
        )

    order = np.array([other_index[tid] for tid in ids_ref], dtype=np.int64)

    return tuple(arr[order] for arr in arrays_other)


def weighted_quantile(values, weights, quantiles):
    """
    Weighted quantile using interpolation on the weighted CDF.

    values, weights: 1D arrays after masking
    quantiles: fractions, e.g. [0.14, 0.50, 0.86]
    """
    values = np.asarray(values, dtype=float)
    weights = np.asarray(weights, dtype=float)
    quantiles = np.asarray(quantiles, dtype=float)

    if values.size == 0:
        return np.full(quantiles.shape, np.nan, dtype=float)

    mask = np.isfinite(values) & np.isfinite(weights) & (weights > 0.0)
    values = values[mask]
    weights = weights[mask]

    if values.size == 0 or np.sum(weights) <= 0.0:
        return np.full(quantiles.shape, np.nan, dtype=float)

    idx = np.argsort(values)
    x = values[idx]
    w = weights[idx]

    # CDF at bin centers. This behaves reasonably for discrete tracer weights.
    cdf = (np.cumsum(w) - 0.5 * w) / np.sum(w)

    return np.interp(quantiles, cdf, x, left=x[0], right=x[-1])


def percentile_label(p):
    """
    Make a compact column label from a percentile value.
    Examples:
      14    -> p14
      50    -> p50
      2.5   -> p2p5
      97.5  -> p97p5
    """
    return "p" + f"{p:g}".replace(".", "p")


def add_range_condition(mask, values, vmin=None, vmax=None, open_upper=False):
    """
    Add finite-value and optional lower/upper cuts.
    """
    cond = np.isfinite(values)
    if vmin is not None:
        cond &= (values >= vmin)
    if vmax is not None:
        if open_upper:
            cond &= (values < vmax)
        else:
            cond &= (values <= vmax)
    return mask & cond

def make_aligned_header(header_cols, widths, comment="# "):
    parts = []
    for i, (name, width) in enumerate(zip(header_cols, widths)):
        if i == 0:
            parts.append(f"{name:>{width - len(comment)}s}")
        else:
            parts.append(f"{name:>{width}s}")
    return comment + " ".join(parts)


def write_percentile_table(output_file, quantity_name, temps, values, weights, base_mask, percentiles, quantiles):
    rows = []

    for j, temp in enumerate(temps):
        value_j = values[:, j]
        mask = base_mask & np.isfinite(value_j)
        qvals = weighted_quantile(value_j[mask], weights[mask], quantiles)
        selected_mass = np.sum(weights[mask])
        n_selected = np.count_nonzero(mask)
        rows.append(np.concatenate(([temp], qvals, [selected_mass, n_selected])))

    rows = np.array(rows)

    q_labels = [f"{quantity_name}_{percentile_label(p)}_mass_weighted" for p in percentiles]
    header_cols = ["T"] + q_labels + ["selected_mass", "n_selected"]

    widths = [22]
    widths += [max(28, len(label) + 2) for label in q_labels]
    widths += [22, 12]

    header = make_aligned_header(header_cols, widths, comment="# ")

    fmt = " ".join(
        [f"%{widths[0]}.8e"]
        + [f"%{widths[i]}.8e" for i in range(1, 1 + len(percentiles))]
        + [f"%{widths[-2]}.8e", f"%{widths[-1]}d"]
    )

    np.savetxt(output_file, rows, header=header, fmt=fmt, comments="")
    print(f"Wrote: {output_file}")

def make_output_name(output_prefix, quantity_name):
    prefix = Path(output_prefix)
    if prefix.suffix in [".txt", ".dat"]:
        prefix = prefix.with_suffix("")
    return str(prefix.parent / f"{prefix.name}_{quantity_name}.txt")

def make_figure_selected(
    ye,
    s,
    vr,
    ye_min,
    ye_max,
    s_min,
    s_max,
    vr_min,
    vr_max,
    prefix,
    highlight_groups=None,
):
    """
    Plot Ye-s, Ye-vr, and s-vr distributions.

    Parameters
    ----------
    highlight_groups : sequence of dict, optional
        Each dictionary must contain

          mask  : boolean array with the same length as ye
          label : legend label

        It may optionally contain ``color``. When color is omitted,
        colors are assigned from Matplotlib's default color cycle.
        The same group uses the same color in all three figures.
    """

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
        'text.usetex'         : True
    }
    plt.rcParams.update(params)

    ye = np.asarray(ye)
    s = np.asarray(s)
    vr = np.asarray(vr)

    if not (ye.shape == s.shape == vr.shape):
        raise ValueError(
            "ye, s, and vr must have identical shapes: "
            f"{ye.shape}, {s.shape}, {vr.shape}"
        )

    if highlight_groups is None:
        highlight_groups = []

    default_colors = plt.rcParams["axes.prop_cycle"].by_key().get(
        "color", [f"C{i}" for i in range(10)]
    )

    prepared_groups = []
    for igroup, group in enumerate(highlight_groups):
        if "mask" not in group:
            raise ValueError(f"highlight_groups[{igroup}] has no 'mask'.")

        group_mask = np.asarray(group["mask"], dtype=bool)
        if group_mask.shape != ye.shape:
            raise ValueError(
                f"highlight_groups[{igroup}]['mask'] has shape "
                f"{group_mask.shape}, expected {ye.shape}."
            )

        prepared_groups.append(
            {
                "mask": group_mask,
                "label": group.get("label", f"group {igroup + 1}"),
                "color": group.get(
                    "color", default_colors[igroup % len(default_colors)]
                ),
            }
        )

    def draw_points(ax, x, y):
        finite = np.isfinite(x) & np.isfinite(y)

        # All tracers are shown faintly in the background.
        ax.scatter(
            x[finite],
            y[finite],
            marker=".",
            s=3,
            color="0.75",
            linewidths=0.0,
            zorder=1,
        )

        # The --inspect-T groups are drawn on top, with one color per bin.
        for group in prepared_groups:
            selected = finite & group["mask"]
            ax.scatter(
                x[selected],
                y[selected],
                marker="o",
                s=10,
                color=group["color"],
                linewidths=0.0,
                alpha=0.85,
                label=group["label"],
                zorder=3,
            )

    def add_legend(ax):
        if prepared_groups:
            ax.legend(
                fontsize=12,
                frameon=False,
                ncol=1,
                loc="upper left",
            )

    ye_min_glo = 0.0
    ye_max_glo = 0.6
    s_min_glo = 1.0
    s_max_glo = 3e2
    vr_min_glo = 1e8
    vr_max_glo = 3e10

    ye0 = ye_min if ye_min is not None else ye_min_glo
    ye1 = ye_max if ye_max is not None else ye_max_glo
    s0 = s_min if s_min is not None else s_min_glo
    s1 = s_max if s_max is not None else s_max_glo
    vr0 = vr_min if vr_min is not None else vr_min_glo
    vr1 = vr_max if vr_max is not None else vr_max_glo

    size_fig = 10.0

    fig_Ye_S = plt.figure(figsize=(size_fig, size_fig*0.625))
    ax_Ye_S = fig_Ye_S.add_subplot(111)

    draw_points(ax_Ye_S, ye, s)

    rect = Rectangle(
        (ye0, s0),
        ye1 - ye0,
        s1 - s0,
        fill=False,
        edgecolor="red",
        linewidth=2,
        zorder=2,
    )
    ax_Ye_S.add_patch(rect)

    fig_Ye_S.subplots_adjust(left=0.15,bottom=0.15,right=0.95, top=0.95)
    add_legend(ax_Ye_S)

    ax_Ye_S.set_xlabel("$Y_\\mathrm{e}$")
    ax_Ye_S.set_xlim(ye_min_glo, ye_max_glo)
    ax_Ye_S.set_xscale("linear")

    ax_Ye_S.set_ylabel("$s/k_\\mathrm{B}$")
    ax_Ye_S.set_ylim(s_min_glo, s_max_glo)
    ax_Ye_S.set_yscale("log")
    ax_Ye_S.yaxis.set_major_locator(ticker.LogLocator(numticks=10))
    ax_Ye_S.yaxis.set_minor_locator(
        ticker.LogLocator(numticks=10, subs=(.1,.2,.3,.4,.5,.6,.7,.8,.9))
    )

    fig_Ye_S.savefig("Ye-S-%s.png" % (prefix))


    fig_Ye_vr = plt.figure(figsize=(size_fig, size_fig*0.625))
    ax_Ye_vr = fig_Ye_vr.add_subplot(111)

    draw_points(ax_Ye_vr, ye, vr)

    rect = Rectangle(
        (ye0, vr0),
        ye1 - ye0,
        vr1 - vr0,
        fill=False,
        edgecolor="red",
        linewidth=2,
        zorder=2,
    )
    ax_Ye_vr.add_patch(rect)

    fig_Ye_vr.subplots_adjust(left=0.15,bottom=0.15,right=0.95, top=0.95)
    add_legend(ax_Ye_vr)

    ax_Ye_vr.set_xlabel("$Y_\\mathrm{e}$")
    ax_Ye_vr.set_xlim(ye_min_glo, ye_max_glo)
    ax_Ye_vr.set_xscale("linear")

    ax_Ye_vr.set_ylabel("$v^r$ (cm\\,s$^{-1}$)")
    ax_Ye_vr.set_ylim(vr_min_glo, vr_max_glo)
    ax_Ye_vr.set_yscale("log")
    ax_Ye_vr.yaxis.set_major_locator(ticker.LogLocator(numticks=10))
    ax_Ye_vr.yaxis.set_minor_locator(
        ticker.LogLocator(numticks=10, subs=(.1,.2,.3,.4,.5,.6,.7,.8,.9))
    )

    fig_Ye_vr.savefig("Ye-vr-%s.png" % (prefix))


    fig_S_vr = plt.figure(figsize=(size_fig, size_fig*0.625))
    ax_S_vr = fig_S_vr.add_subplot(111)

    draw_points(ax_S_vr, s, vr)

    rect = Rectangle(
        (s0, vr0),
        s1 - s0,
        vr1 - vr0,
        fill=False,
        edgecolor="red",
        linewidth=2,
        zorder=2,
    )
    ax_S_vr.add_patch(rect)

    fig_S_vr.subplots_adjust(left=0.15,bottom=0.15,right=0.95, top=0.95)
    add_legend(ax_S_vr)

    ax_S_vr.set_xlabel("$s/k_\\mathrm{B}$")
    ax_S_vr.set_xlim(s_min_glo, s_max_glo)
    ax_S_vr.set_xscale("log")
    ax_S_vr.xaxis.set_major_locator(ticker.LogLocator(numticks=15))
    ax_S_vr.xaxis.set_minor_locator(
        ticker.LogLocator(numticks=10, subs=(.1,.2,.3,.4,.5,.6,.7,.8,.9))
    )

    ax_S_vr.set_ylabel("$v^r$ (cm\\,s$^{-1}$)")
    ax_S_vr.set_ylim(vr_min_glo, vr_max_glo)
    ax_S_vr.set_yscale("log")
    ax_S_vr.yaxis.set_major_locator(ticker.LogLocator(numticks=10))
    ax_S_vr.yaxis.set_minor_locator(
        ticker.LogLocator(numticks=10, subs=(.1,.2,.3,.4,.5,.6,.7,.8,.9))
    )

    fig_S_vr.savefig("S-vr-%s.png" % (prefix))

def main():
    parser = argparse.ArgumentParser(
        description=(
            "Compute mass-weighted percentiles of a temperature-dependent tracer quantity "
            "using only tracers selected by final Ye, entropy, and radial velocity from stat_traj.dat."
        )
    )

    parser.add_argument("--time_file", type=str, required=False, default="traj_time", help="Temperature table to percentile, e.g. traj_texp.txt")
    parser.add_argument("--texp_file", type=str, required=False, default="traj_texp", help="Temperature table to percentile, e.g. traj_texp.txt")
    parser.add_argument("--rcap_file", type=str, required=False, default="traj_caprate", help="Temperature table to percentile, e.g. traj_texp.txt")
    parser.add_argument("--ye_file", type=str, required=False, default="traj_Ye", help="Temperature table to percentile, e.g. traj_texp.txt")
    parser.add_argument("--yecap_file", type=str, required=False, default="traj_yecap", help="Temperature table to percentile, e.g. traj_texp.txt")
    parser.add_argument("--eta_file", type=str, required=False, default="traj_eta", help="Temperature table to percentile, e.g. traj_texp.txt")
    parser.add_argument("--r_file", type=str, required=False, default="traj_r", help="Temperature table to percentile, e.g. traj_texp.txt")
    parser.add_argument("--vr_file", type=str, required=False, default="traj_vel", help="Temperature table to percentile, e.g. traj_texp.txt")
    parser.add_argument("--avtexp_file", type=str, required=False, default="traj_avtexp", help="Temperature table to percentile, e.g. traj_texp.txt")
    parser.add_argument("--stat_file", type=str, required=False, default="stat_traj.dat", help="stat_traj.dat")

    parser.add_argument(
        "--ye-min", type=float, default=None,
        help="Minimum Ye_fin. Omit for no lower cut."
    )
    parser.add_argument(
        "--ye-max", type=float, default=None,
        help="Maximum Ye_fin. Omit for no upper cut."
    )
    parser.add_argument(
        "--s-min", type=float, default=None,
        help="Minimum s_fin [k_b/nuc]. Omit for no lower cut."
    )
    parser.add_argument(
        "--s-max", type=float, default=None,
        help="Maximum s_fin [k_b/nuc]. Omit for no upper cut."
    )
    parser.add_argument(
        "--vr-min", type=float, default=None,
        help="Minimum v^r_fin [cm/s]. Omit for no lower cut."
    )
    parser.add_argument(
        "--vr-max", type=float, default=None,
        help="Maximum v^r_fin [cm/s]. Omit for no upper cut."
    )

    parser.add_argument(
        "--open-upper",
        action="store_true",
        help="Use < for all upper bounds instead of <=."
    )
    parser.add_argument(
        "--percentiles",
        type=float,
        nargs="+",
        default=[15.865, 50.0, 84.135], # 1sigma
        help=(
            "Percentiles to compute in percent, e.g. "
            "--percentiles 5 50 95 or --percentiles 2.5 50 97.5. "
            "Default: 14 50 86"
        ),
    )
    parser.add_argument(
        "--weight-source",
        choices=["stat", "value"],
        default="stat",
        help="Use masses from stat_traj.dat or from the value table as weights. Default: stat."
    )
    # parser.add_argument(
    #     "-o", "--output",
    #     default="percentiles_selected_by_final_quantities.txt",
    #     help="Output filename"
    # )
    parser.add_argument(
        "-o", "--output-prefix",
        default="percentiles_selected_by_final_quantities",
        help=(
            "Output prefix. The script writes "
            "<prefix>_texp.txt, <prefix>_ye.txt, <prefix>_yecap.txt, and <prefix>_caprate.txt. "
            "If a .txt/.dat suffix is given, it is stripped."
        ),
    )


    args = parser.parse_args()

    percentiles = np.array(args.percentiles, dtype=float)
    if np.any((percentiles < 0.0) | (percentiles > 100.0)):
        raise ValueError("All percentiles must be between 0 and 100.")
    quantiles = percentiles / 100.0

    temps, ids, mass, time = read_temperature_table(args.time_file)
    _, _, _, texp = read_temperature_table(args.texp_file)
    _, _, _, caprate = read_temperature_table(args.rcap_file)
    _, _, _, ye = read_temperature_table(args.ye_file)
    _, _, _, yecap = read_temperature_table(args.yecap_file)
    _, _, _, eta = read_temperature_table(args.eta_file)
    _, _, _, r = read_temperature_table(args.r_file)
    _, _, _, vr= read_temperature_table(args.vr_file)
    _, _, _, avtexp = read_temperature_table(args.avtexp_file)
    (
        ids_stat,
        mass_stat_raw,
        ye_fin_raw,
        s_fin_raw,
        vr_fin_raw,
        ye_5GK_raw,
        s_5GK_raw,
        bernoulli_raw,
        r_fin_raw,
    ) = read_stat_traj(args.stat_file)

    (
        mass_stat,
        ye_fin,
        s_fin,
        vr_fin,
        ye_5GK,
        s_5GK,
        bernoulli,
        r_fin,
    ) = align_to_reference_ids(
        ids,
        ids_stat,
        mass_stat_raw,
        ye_fin_raw,
        s_fin_raw,
        vr_fin_raw,
        ye_5GK_raw,
        s_5GK_raw,
        bernoulli_raw,
        r_fin_raw,
        label=args.stat_file,
    )

    if not np.allclose(mass, mass_stat, rtol=1.0e-10, atol=0.0):
        max_rel = np.nanmax(
            np.abs(mass - mass_stat) / np.maximum(np.abs(mass_stat), 1.0)
        )
        print(
            "Warning: masses in value file and stat file are not identical. "
            f"max relative difference = {max_rel:.3e}"
        )

    if args.weight_source == "stat":
        weights = mass_stat
    else:
        weights = mass

    base_mask = np.isfinite(weights) & (weights > 0.0)
    
    base_mask = add_range_condition(
        base_mask, ye_5GK, args.ye_min, args.ye_max, open_upper=args.open_upper
    )
    base_mask = add_range_condition(
        base_mask, s_5GK, args.s_min, args.s_max, open_upper=args.open_upper
    )
    base_mask = add_range_condition(
        base_mask, vr_fin, args.vr_min, args.vr_max, open_upper=args.open_upper
    )

    fn_out = args.output_prefix + "_ids.txt"
    with open(fn_out, mode="w") as f:
        f.write("# IDs for selected trajectories\n")
        for ip in ids[base_mask]:
            f.write("%d\n" % (ip))
            pass
        pass
    make_figure_selected(ye_5GK, s_5GK, vr_fin, 
                         args.ye_min, args.ye_max,
                         args.s_min, args.s_max,
                         args.vr_min, args.vr_max,
                         args.output_prefix)
    
    n_base = np.count_nonzero(base_mask)
    mass_base = np.sum(weights[base_mask])
    
    quantities = {
        "ye": ye,
        "texp": texp,
        "yecap": yecap,
        "caprate": caprate,
        "eta": eta,
        "r": r,
        "vr": vr,
        "time": time,
        "avtexp": avtexp,
    }
    
    for quantity_name, values in quantities.items():
        output_file = make_output_name(args.output_prefix, quantity_name)
        write_percentile_table(
            output_file=output_file,
            quantity_name=quantity_name,
            temps=temps,
            values=values,
            weights=weights,
            base_mask=base_mask,
            percentiles=percentiles,
            quantiles=quantiles,
        )
        
    # rows = []

    # for j, temp in enumerate(temps):
    #     texp_j = texp[:, j]
    #     ye_j = ye[:, j]
    #     caprate_j = caprate[:, j]
    #     yecap_j = yecap[:, j]

    #     # The final-quantity selection is fixed for all temperatures.
    #     # The finite-value mask can depend on temperature.
    #     mask = base_mask & np.isfinite(ye_j)

    #     qvals = weighted_quantile(texp_j[mask], weights[mask], quantiles)
    #     selected_mass = np.sum(weights[mask])
    #     n_selected = np.count_nonzero(mask)

    #     rows.append(np.concatenate(([temp], qvals, [selected_mass, n_selected])))

    # rows = np.array(rows)

    # q_labels = [f"value_{percentile_label(p)}_mass_weighted" for p in percentiles]
    # header_cols = ["T"] + q_labels + ["selected_mass", "n_selected"]
    
    # widths = [22] + [28] * len(percentiles) + [22, 12]
    
    # # header = "# " + " ".join(
    # #     f"{name:>{w}s}" for name, w in zip(header_cols, widths)
    # # )

    # header_parts = []
    
    # for i, (name, w) in enumerate(zip(header_cols, widths)):
    #     if i == 0:
    #         # "# " の2文字分だけ、最初の列幅を短くする
    #         header_parts.append(f"{name:>{w-2}s}")
    #     else:
    #         header_parts.append(f"{name:>{w}s}")
    #         pass
    #     pass
            
    # header = "# " + " ".join(header_parts)
    
    # fmt = " ".join(
    #     [f"%{widths[0]}.8e"]
    #     + [f"%{w}.8e" for w in widths[1:-2]]
    #     + [f"%{widths[-2]}.8e", f"%{widths[-1]}d"]
    # )
    
    # np.savetxt(
    #     args.output,
    #     rows,
    #     header=header,
    #     fmt=fmt,
    #     comments="",
    # )

    # # q_labels = [f"value_{percentile_label(p)}_mass_weighted" for p in percentiles]
    # # header_cols = ["T"] + q_labels + ["selected_mass", "n_selected"]
    
    # # fmt = ["%.8e"] * (1 + len(percentiles) + 1) + ["%d"]

    # # np.savetxt(
    # #     args.output,
    # #     rows,
    # #     header=" ".join(header_cols),
    # #     fmt=fmt,
    # # )

    # print(f"Wrote: {args.output}")

    print(f"Base selected tracers: n = {n_base}, mass = {mass_base:.8e} g")
    print("Cuts:")
    print(f"  Ye_fin: {args.ye_min} <= Ye_fin {'<' if args.open_upper else '<='} {args.ye_max}")
    print(f"  s_fin : {args.s_min} <= s_fin  {'<' if args.open_upper else '<='} {args.s_max}")
    print(f"  vr_fin: {args.vr_min} <= vr_fin {'<' if args.open_upper else '<='} {args.vr_max} cm/s")
    print("Percentiles:", " ".join(f"{p:g}" for p in percentiles))


if __name__ == "__main__":
    main()
    plt.show()
