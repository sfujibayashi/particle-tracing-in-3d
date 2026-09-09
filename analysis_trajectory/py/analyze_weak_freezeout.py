import argparse
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import matplotlib.font_manager as fm
import matplotlib.cm as cm
import sys
import os
import re
import shlex

import tracer_analysis_utils as stat


MeV_to_K=1.16e10

params={
    'font.size'           : 24.0     ,
#    'font.family'         : 'DejaVu Sans'     ,
    'font.family'         : 'Times New Roman'  ,
    'xtick.major.size'    : 2        ,
    'xtick.major.width'   : 1.5      ,
    'xtick.labelsize'     : 24.0     ,
    'xtick.direction'     : 'in'     ,
    'ytick.major.size'    : 2        ,
    'ytick.major.width'   : 1.5      ,
    'ytick.labelsize'     : 24.0     ,
    'ytick.direction'     : 'in'     ,
    'xtick.major.pad'     : 2        ,
    'xtick.minor.pad'     : 2        ,
    'ytick.major.pad'     : 2        ,
    'ytick.minor.pad'     : 2        ,
#    'savefig.dpi'         : 601      ,
    'axes.linewidth'      : 1.5      ,
    'text.usetex'         : True    }

plt.rcParams.update(params)

# Integrated here so this workflow does not depend on a separate
# make_figure_selected_5GK.py file.
def make_figure_selected_5GK(
    stat_file,
    prefix,
    highlight_groups=None,
    ids_ref=None,
    background_mask=None,
    output_ext="png",
    ids_ref_label="ids_ref",
    ye_range=None,
    s_range=None,
    texp_range=None,
    annotation_text=None,
    annotation_fontsize=11,
    annotation_loc="upper left",
):
    """
    Read s_5GK, Ye_5GK, and t_exp_5GK from stat_traj.dat and make

      1. Ye_5GK    vs s_5GK
      2. Ye_5GK    vs t_exp_5GK
      3. s_5GK     vs t_exp_5GK

    scatter plots with three distinct layers:

      1. all rows in stat_traj.dat:
         small gray points;

      2. rows whose PIDs are listed in ids_ref:
         open black circles, optionally restricted by background_mask;

      3. highlight_groups:
         larger filled colored circles with black edges.

    Parameters
    ----------
    stat_file : str or pathlib.Path
        stat_traj.dat filename.

    prefix : str or pathlib.Path
        Output prefix. For example, prefix="hist_freezeout" produces

          Ye-S_5GK-hist_freezeout.png
          Ye-texp_5GK-hist_freezeout.png
          S-texp_5GK-hist_freezeout.png

    highlight_groups : sequence of dict, optional
        Tracer groups to overplot. Each dictionary must contain

          "mask"  : boolean array
          "label" : legend label

        and may optionally contain

          "color" : Matplotlib color

        The same group is drawn with the same color in all three figures.

    ids_ref : array-like of int, optional
        Reference PID order. When supplied, quantities read from stat_traj.dat
        are reordered to match ids_ref for the second and third layers.

        Irrespective of ids_ref, every finite row of stat_traj.dat is shown
        as a small gray point in the first layer.

    background_mask : array-like of bool, optional
        Mask selecting which ids_ref tracers are shown in the second layer.
        It must be in ids_ref order when ids_ref is supplied. The default
        shows all ids_ref tracers.

        This mask also restricts highlight_groups, preserving the previous
        behavior that highlighted points are a subset of the selected sample.

    output_ext : str, optional
        Output extension, without a leading period. Default: "png".

    ids_ref_label : str or None, optional
        Legend label for the ids_ref layer. Set to None to omit that layer
        from the legend. Default: "ids_ref".

    annotation_text : str or None, optional
        Multiline text drawn inside all three scatter plots. Use this to
        record the cuts that produced background_mask. Default: None.

    annotation_fontsize : float, optional
        Font size of annotation_text. Default: 11.

    annotation_loc : {"upper left", "upper right", "lower left", "lower right"}
        Position of annotation_text in axes coordinates.

    Returns
    -------
    dict
        Reordered arrays and output filenames.
    """
    from pathlib import Path
    import numpy as np
    import matplotlib.pyplot as plt
    import matplotlib.ticker as ticker

    # stat_traj.dat, 0-based column indices:
    #   0  : id
    #   13 : s_5GK [k_b/nuc]
    #   14 : Ye_5GK
    #   15 : t_exp_5GK [s]
    arr = np.loadtxt(
        stat_file,
        comments="#",
        usecols=(0, 13, 14, 15),
        # usecols=(0, 10, 11, 15),
        ndmin=2,
    )

    ids_stat = arr[:, 0].astype(np.int64)
    s_5GK_stat = arr[:, 1]
    ye_5GK_stat = arr[:, 2]
    t_exp_5GK_stat = arr[:, 3]

    if np.unique(ids_stat).size != ids_stat.size:
        raise ValueError(f"{stat_file}: duplicate tracer IDs are present.")

    ids_ref_supplied = ids_ref is not None

    if ids_ref_supplied:
        ids_ref = np.asarray(ids_ref, dtype=np.int64)

        if ids_ref.ndim != 1:
            raise ValueError(
                f"ids_ref must be one-dimensional; got shape {ids_ref.shape}."
            )
        if np.unique(ids_ref).size != ids_ref.size:
            raise ValueError("ids_ref contains duplicate tracer IDs.")

        stat_index = {pid: i for i, pid in enumerate(ids_stat)}
        missing = [pid for pid in ids_ref if pid not in stat_index]
        if missing:
            raise ValueError(
                f"{len(missing)} tracer IDs are missing from {stat_file}. "
                f"First missing PID: {missing[0]}"
            )

        order = np.fromiter(
            (stat_index[pid] for pid in ids_ref),
            dtype=np.int64,
            count=ids_ref.size,
        )

        ids_plot = ids_stat[order]
        s_5GK = s_5GK_stat[order]
        ye_5GK = ye_5GK_stat[order]
        t_exp_5GK = t_exp_5GK_stat[order]
    else:
        ids_plot = ids_stat.copy()
        s_5GK = s_5GK_stat.copy()
        ye_5GK = ye_5GK_stat.copy()
        t_exp_5GK = t_exp_5GK_stat.copy()

    ntracer = ids_plot.size

    background_mask_supplied = background_mask is not None

    if background_mask is None:
        background_mask = np.ones(ntracer, dtype=bool)
    else:
        background_mask = np.asarray(background_mask, dtype=bool)
        if background_mask.shape != (ntracer,):
            raise ValueError(
                "background_mask has shape "
                f"{background_mask.shape}; expected {(ntracer,)}."
            )

    # When neither ids_ref nor background_mask is supplied, avoid drawing
    # every stat_traj.dat point twice.
    draw_ids_ref_layer = ids_ref_supplied or background_mask_supplied

    if highlight_groups is None:
        highlight_groups = []

    default_colors = plt.rcParams["axes.prop_cycle"].by_key().get(
        "color", [f"C{i}" for i in range(10)]
    )

    groups = []
    for igroup, group in enumerate(highlight_groups):
        if "mask" not in group:
            raise ValueError(f"highlight_groups[{igroup}] has no 'mask'.")

        group_mask = np.asarray(group["mask"], dtype=bool)
        if group_mask.shape != (ntracer,):
            raise ValueError(
                f"highlight_groups[{igroup}]['mask'] has shape "
                f"{group_mask.shape}; expected {(ntracer,)}."
            )

        groups.append(
            {
                "mask": group_mask,
                "label": group.get("label", f"group {igroup + 1}"),
                "color": group.get(
                    "color",
                    default_colors[igroup % len(default_colors)],
                ),
                "marker": group.get("marker", "."),
                "size": float(group.get("size", 5.0)),
            }
        )

    prefix_path = Path(prefix)
    output_dir = prefix_path.parent
    output_dir.mkdir(parents=True, exist_ok=True)
    stem = prefix_path.name
    ext = output_ext.lstrip(".")

    output_files = {
        "Ye-S": output_dir / f"Ye-S_5GK-{stem}.{ext}",
        "Ye-texp": output_dir / f"Ye-texp_5GK-{stem}.{ext}",
        "S-texp": output_dir / f"S-texp_5GK-{stem}.{ext}",
    }

    def finite_for_scale(x, y, xscale, yscale):
        finite = np.isfinite(x) & np.isfinite(y)
        if xscale == "log":
            finite &= x > 0.0
        if yscale == "log":
            finite &= y > 0.0
        return finite

    def draw_annotation(ax):
        if not annotation_text:
            return

        locations = {
            "upper left": (0.02, 0.98, "left", "top"),
            "upper right": (0.98, 0.98, "right", "top"),
            "lower left": (0.02, 0.02, "left", "bottom"),
            "lower right": (0.98, 0.02, "right", "bottom"),
        }
        if annotation_loc not in locations:
            raise ValueError(
                "annotation_loc must be one of "
                f"{tuple(locations)}; got {annotation_loc!r}."
            )

        x, y, ha, va = locations[annotation_loc]
        ax.text(
            x,
            y,
            annotation_text,
            transform=ax.transAxes,
            ha=ha,
            va=va,
            fontsize=annotation_fontsize,
            zorder=20,
            bbox=dict(
                boxstyle="round,pad=0.35",
                facecolor="white",
                edgecolor="0.7",
                alpha=0.85,
            ),
        )

    def plot_scatter(
        x_stat,
        y_stat,
        x_ref,
        y_ref,
        xlabel,
        ylabel,
        output_file,
        xscale="linear",
        yscale="linear",
        xlim=None,
        ylim=None,
    ):
        fig = plt.figure(figsize=(10.0, 6.25))
        ax = fig.add_subplot(111)

        # Layer 1: all finite rows of stat_traj.dat.
        finite_stat = finite_for_scale(
            x_stat,
            y_stat,
            xscale,
            yscale,
        )
        ax.scatter(
            x_stat[finite_stat],
            y_stat[finite_stat],
            marker=".",
            s=3,
            color="0.75",
            linewidths=0.0,
            alpha=0.8,
            zorder=1,
        )

        finite_ref = finite_for_scale(
            x_ref,
            y_ref,
            xscale,
            yscale,
        )

        # Layer 2: ids_ref, optionally restricted by background_mask.
        selected_reference = finite_ref & background_mask
        if draw_ids_ref_layer and np.any(selected_reference):
            ax.scatter(
                x_ref[selected_reference],
                y_ref[selected_reference],
                marker=".",
                s=3,
                facecolors="none",
                edgecolors="r",
                linewidths=0.7,
                alpha=0.9,
                label=ids_ref_label,
                zorder=2,
            )

        # Layer 3: highlight_groups, drawn above both previous layers.
        for group in groups:
            selected = (
                finite_ref
                & background_mask
                & group["mask"]
            )
            ax.scatter(
                x_ref[selected],
                y_ref[selected],
                marker=group["marker"],
                s=group["size"],
                facecolors=group["color"],
                #edgecolors="black",
                linewidths=0.6,
                alpha=0.95,
                label=group["label"],
                zorder=4,
            )

        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.set_xscale(xscale)
        ax.set_yscale(yscale)

        if xlim is not None:
            ax.set_xlim(xlim)
        
        if ylim is not None:
            ax.set_ylim(ylim)

        ax.tick_params(
            axis="both",
            which="major",
            direction="out",
            width=1.5,
            length=6,
        )
        ax.tick_params(
            axis="both",
            which="minor",
            direction="out",
            width=1.5,
            length=3,
        )

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

        draw_annotation(ax)

        handles, labels = ax.get_legend_handles_labels()
        if handles:
            ax.legend(
                fontsize=12,
                frameon=False,
                ncol=1,
                loc="best",
            )

        fig.subplots_adjust(
            left=0.15,
            bottom=0.15,
            right=0.95,
            top=0.95,
        )
        fig.savefig(output_file)
        #plt.close(fig)

    plot_scatter(
        ye_5GK_stat,
        s_5GK_stat,
        ye_5GK,
        s_5GK,
        r"$Y_{\rm e}(5\,{\rm GK})$",
        r"$s(5\,{\rm GK})/k_{\rm B}$",
        output_files["Ye-S"],
        xscale="linear",
        yscale="log",
        xlim=ye_range,
        ylim=s_range,
    )

    plot_scatter(
        ye_5GK_stat,
        t_exp_5GK_stat,
        ye_5GK,
        t_exp_5GK,
        r"$Y_{\rm e}(5\,{\rm GK})$",
        r"$t_{\rm exp}(5\,{\rm GK})\ [{\rm s}]$",
        output_files["Ye-texp"],
        xscale="linear",
        yscale="log",
        xlim=ye_range,
        ylim=texp_range,
    )

    plot_scatter(
        s_5GK_stat,
        t_exp_5GK_stat,
        s_5GK,
        t_exp_5GK,
        r"$s(5\,{\rm GK})/k_{\rm B}$",
        r"$t_{\rm exp}(5\,{\rm GK})\ [{\rm s}]$",
        output_files["S-texp"],
        xscale="log",
        yscale="log",
        xlim=s_range,
        ylim=texp_range,
    )

    return {
        "ids_stat": ids_stat,
        "s_5GK_stat": s_5GK_stat,
        "Ye_5GK_stat": ye_5GK_stat,
        "t_exp_5GK_stat": t_exp_5GK_stat,
        "ids": ids_plot,
        "s_5GK": s_5GK,
        "Ye_5GK": ye_5GK,
        "t_exp_5GK": t_exp_5GK,
        "background_mask": background_mask,
        "output_files": {
            key: str(value) for key, value in output_files.items()
        },
    }

parser = argparse.ArgumentParser()
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
    "--r-fin-min",
    type=float,
    default=None,
    help="Minimum final coordinate radius [cm]. Default: 1e8.",
)
parser.add_argument(
    "--r-fin-max",
    type=float,
    default=None,
    help="Minimum final coordinate radius [cm]. Default: 1e8.",
)
parser.add_argument(
    "--time_fin_min",
    "--time-fin-min",
    dest="time_fin_min",
    type=float,
    default=None,
    help=(
        "Minimum time_fin [s]. Omit for no lower cut. "
        "The underscore and hyphen spellings are equivalent."
    ),
)
parser.add_argument(
    "--time_fin_max",
    "--time-fin-max",
    dest="time_fin_max",
    type=float,
    default=None,
    help=(
        "Maximum time_fin [s]. Omit for no upper cut. "
        "By default the upper bound is inclusive; use --open-upper "
        "to make it strict."
    ),
)
parser.add_argument(
    "--time-fin-column",
    type=int,
    default=None,
    help=(
        "Optional 1-based column number of time_fin in stat_traj.dat. "
        "Normally it is detected from the header automatically."
    ),
)


parser.add_argument(
    "--time_rfl_min",
    "--time-rfl-min",
    dest="time_rfl_min",
    type=float,
    default=0.0,
    help=(
        "Minimum time_rfl [s]. Omit for no lower cut. "
        "The underscore and hyphen spellings are equivalent."
    ),
)
parser.add_argument(
    "--time_rfl_max",
    "--time-rfl-max",
    dest="time_rfl_max",
    type=float,
    default=None,
    help=(
        "Maximum time_rfl [s]. Omit for no upper cut. "
        "By default the upper bound is inclusive; use --open-upper "
        "to make it strict."
    ),
)
parser.add_argument(
    "--time-rfl-column",
    type=int,
    default=None,
    help=(
        "Optional 1-based column number of time_rfl in stat_traj.dat. "
        "Normally it is detected from the header automatically."
    ),
)

parser.add_argument(
    "--theta_fin_min",
    "--theta-fin-min",
    "--theta-min",
    dest="theta_fin_min",
    type=float,
    default=None,
    help=(
        "Minimum final-position angle [deg]. Omit for no lower cut. "
        "The angle definition is selected with --theta-mode."
    ),
)
parser.add_argument(
    "--theta_fin_max",
    "--theta-fin-max",
    "--theta-max",
    dest="theta_fin_max",
    type=float,
    default=None,
    help=(
        "Maximum final-position angle [deg]. Omit for no upper cut. "
        "By default the upper bound is inclusive; use --open-upper "
        "to make it strict."
    ),
)
parser.add_argument(
    "--theta-mode",
    choices=("axis", "polar", "latitude"),
    default="axis",
    help=(
        "Definition of the angle made from x_fin, y_fin, z_fin. "
        "'axis' (default): angle from the nearest +/-z axis, "
        "atan2(sqrt(x^2+y^2), abs(z)), range 0--90 deg. "
        "'polar': standard polar angle from +z, "
        "atan2(sqrt(x^2+y^2), z), range 0--180 deg. "
        "'latitude': signed angle from the equatorial plane, "
        "atan2(z, sqrt(x^2+y^2)), range -90--90 deg."
    ),
)

parser.add_argument(
    "--cut-label-loc",
    choices=("upper left", "upper right", "lower left", "lower right"),
    default="upper left",
    help=(
        "Position of the active-cut annotation in all histogram and "
        "5-GK scatter plots. Default: upper left."
    ),
)
parser.add_argument(
    "--cut-label-fontsize",
    type=float,
    default=11.0,
    help="Font size of the active-cut annotation. Default: 11.",
)

parser.add_argument(
    "-o", "--output-prefix",
    default="weakfo",
    help=(
        "Output prefix. The script writes "
        "<prefix>_texp.txt, <prefix>_ye.txt, <prefix>_yecap.txt, and <prefix>_caprate.txt. "
        "If a .txt/.dat suffix is given, it is stripped."
    ),
)

parser.add_argument(
    "--inspect-T",
    type=float,
    nargs="+",
    default=None,
    help=(
        "Inspect the histogram bin containing each specified temperature "
        "[MeV]. When omitted, no T_FO-selected group is written or "
        "highlighted. Do not combine with --inspect-T-range."
    ),
)
parser.add_argument(
    "--inspect-T-range",
    dest="inspect_T_ranges",
    type=float,
    nargs=2,
    action="append",
    metavar=("T_MIN", "T_MAX"),
    default=None,
    help=(
        "Inspect an exact freeze-out-temperature range [T_MIN, T_MAX] MeV. "
        "With --open-upper, the upper boundary is excluded. Repeat this "
        "option to inspect multiple ranges."
    ),
)
parser.add_argument(
    "--inspect-source", choices=("time", "dye", "both"), default="dye",
    help=(
        "Histogram whose particle IDs are written: "
        "'time' for t_weak > t_exp, 'dye' for Delta Ye < 0.01, "
        "or 'both'. Default: dye."
    ),
)


args = parser.parse_args()

if not np.isfinite(args.cut_label_fontsize) or args.cut_label_fontsize <= 0.0:
    parser.error("--cut-label-fontsize must be a positive finite number.")

for option_name, option_value in (
    ("--time_fin_min", args.time_fin_min),
    ("--time_fin_max", args.time_fin_max),
):
    if option_value is not None and not np.isfinite(option_value):
        parser.error(f"{option_name} must be finite.")

if (
    args.time_fin_min is not None
    and args.time_fin_max is not None
    and args.time_fin_min >= args.time_fin_max
):
    parser.error("--time_fin_min must be smaller than --time_fin_max.")

if args.time_fin_column is not None and args.time_fin_column <= 0:
    parser.error("--time-fin-column must be a positive 1-based column number.")

for option_name, option_value in (
    ("--theta-fin-min", args.theta_fin_min),
    ("--theta-fin-max", args.theta_fin_max),
):
    if option_value is not None and not np.isfinite(option_value):
        parser.error(f"{option_name} must be finite.")

theta_bounds = {
    "axis": (0.0, 90.0),
    "polar": (0.0, 180.0),
    "latitude": (-90.0, 90.0),
}
theta_allowed_min, theta_allowed_max = theta_bounds[args.theta_mode]
for option_name, option_value in (
    ("--theta-fin-min", args.theta_fin_min),
    ("--theta-fin-max", args.theta_fin_max),
):
    if option_value is not None and not (
        theta_allowed_min <= option_value <= theta_allowed_max
    ):
        parser.error(
            f"{option_name}={option_value:g} is outside the allowed "
            f"range [{theta_allowed_min:g}, {theta_allowed_max:g}] deg "
            f"for --theta-mode {args.theta_mode}."
        )

if (
    args.theta_fin_min is not None
    and args.theta_fin_max is not None
    and args.theta_fin_min >= args.theta_fin_max
):
    parser.error("--theta-fin-min must be smaller than --theta-fin-max.")

# Omitted temperature-selection options mean that no T_FO-selected
# particle group is constructed. Empty lists also keep downstream loops safe.
if args.inspect_T is None:
    args.inspect_T = []
if args.inspect_T_ranges is None:
    args.inspect_T_ranges = []

if args.inspect_T and args.inspect_T_ranges:
    parser.error(
        "--inspect-T and --inspect-T-range cannot be used together."
    )

if not args.inspect_T and not args.inspect_T_ranges:
    print(
        "No --inspect-T or --inspect-T-range was specified: "
        "T_FO-selected highlighting is disabled, but a particles table "
        "for all tracers passing the active cuts will be written."
    )


def normalize_header_token(token):
    """Normalize a stat_traj.dat header token for column-name matching."""
    token = token.lower().strip()
    # Remove an attached unit, e.g. time_fin[s] -> time_fin.
    token = re.sub(r"\[[^\]]*\]$", "", token)
    token = re.sub(r"\([^)]*\)$", "", token)
    return re.sub(r"[^a-z0-9_]", "", token)


def detect_named_column_from_header(
    filename,
    candidate_names,
):
    """
    Detect a 0-based data-column index from comment headers.

    Supported examples include
      # id mass time_fin x_fin ...
      # 1:id 2:mass 3:time_fin ...
      # 1 id 2 mass 3 time_fin ...

    Returns None when no reliable match is found.
    """
    candidate_names = {
        normalize_header_token(name) for name in candidate_names
    }

    comment_lines = []
    first_data_tokens = None

    with open(filename, "r") as stream:
        for raw_line in stream:
            stripped = raw_line.strip()
            if not stripped:
                continue
            if stripped.startswith("#"):
                comment_lines.append(stripped.lstrip("#").strip())
                continue

            first_data_tokens = stripped.split()
            break

    if first_data_tokens is None:
        raise ValueError(f"{filename}: no data rows were found.")

    n_columns = len(first_data_tokens)

    # Case 1: one header token per data column. Standalone unit tokens such
    # as "[g]" or "[s]" are ignored.
    for line in reversed(comment_lines):
        tokens = line.replace(",", " ").split()
        tokens_without_units = [
            token for token in tokens
            if not re.fullmatch(r"\[[^\]]*\]|\([^)]*\)", token)
        ]
        normalized = [
            normalize_header_token(token) for token in tokens_without_units
        ]

        if len(normalized) == n_columns:
            for index, token in enumerate(normalized):
                if token in candidate_names:
                    return index

    # Case 2: numbered labels such as 3:time_fin or 3 time_fin.
    joined_candidates = "|".join(
        re.escape(name) for name in sorted(candidate_names, key=len, reverse=True)
    )
    numbered_patterns = [
        re.compile(
            rf"(?<!\d)(\d+)\s*[:=.)_-]*\s*({joined_candidates})(?![a-z0-9_])",
            re.IGNORECASE,
        ),
        re.compile(
            rf"({joined_candidates})\s*[:=.(\[]*\s*(\d+)(?!\d)",
            re.IGNORECASE,
        ),
    ]

    for line in reversed(comment_lines):
        normalized_line = line.lower()
        for ipattern, pattern in enumerate(numbered_patterns):
            match = pattern.search(normalized_line)
            if match is None:
                continue

            number_text = match.group(1 if ipattern == 0 else 2)
            column_number = int(number_text)
            if 1 <= column_number <= n_columns:
                return column_number - 1

    return None


def read_time_fin_from_stat(
    filename,
    column_override=None,
):
    """
    Read PID and time_fin from stat_traj.dat.

    column_override is a 1-based data-column number.  When omitted, the
    function searches the comment header for 'time_fin' or 't_fin'.
    """
    if column_override is None:
        time_column = detect_named_column_from_header(
            filename,
            candidate_names=("time_fin", "t_fin", "timefinal"),
        )
        if time_column is None:
            raise ValueError(
                f"Could not identify the time_fin column in {filename}. "
                "Add a header containing the name 'time_fin', or specify "
                "--time-fin-column N using a 1-based column number."
            )
    else:
        time_column = column_override - 1

    try:
        arr = np.loadtxt(
            filename,
            comments="#",
            usecols=(0, time_column),
            ndmin=2,
        )
    except (IndexError, ValueError) as exc:
        raise ValueError(
            f"{filename}: could not read PID and time_fin from 1-based "
            f"column {time_column + 1}."
        ) from exc

    ids_time_fin = arr[:, 0].astype(np.int64)
    time_fin_values = arr[:, 1].astype(float)

    return ids_time_fin, time_fin_values, time_column


def read_time_rfl_from_stat(
    filename,
    column_override=None,
):
    """
    Read PID and time_rfl from stat_traj.dat.

    column_override is a 1-based data-column number.  When omitted, the
    function searches the comment header for 't_last_rfl'.
    """
    if column_override is None:
        time_column = detect_named_column_from_header(
            filename,
            candidate_names=("t_last_rfl"),
        )
        if time_column is None:
            raise ValueError(
                f"Could not identify the time_rfl column in {filename}. "
                "Add a header containing the name 't_last_rfl', or specify "
                "--time-rfl-column N using a 1-based column number."
            )
    else:
        time_column = column_override - 1

    try:
        arr = np.loadtxt(
            filename,
            comments="#",
            usecols=(0, time_column),
            ndmin=2,
        )
    except (IndexError, ValueError) as exc:
        raise ValueError(
            f"{filename}: could not read PID and time_rfl from 1-based "
            f"column {time_column + 1}."
        ) from exc

    ids_time_rfl = arr[:, 0].astype(np.int64)
    time_rfl_values = arr[:, 1].astype(float)

    return ids_time_rfl, time_rfl_values, time_column


def read_final_positions_from_stat(filename):
    """Read PID and final Cartesian coordinates from stat_traj.dat."""
    coordinate_columns = []
    for coordinate_name in ("x_fin", "y_fin", "z_fin"):
        column = detect_named_column_from_header(
            filename,
            candidate_names=(coordinate_name,),
        )
        if column is None:
            raise ValueError(
                f"Could not identify the {coordinate_name} column in "
                f"{filename}. Add a header containing that column name."
            )
        coordinate_columns.append(column)

    try:
        arr = np.loadtxt(
            filename,
            comments="#",
            usecols=(0, *coordinate_columns),
            ndmin=2,
        )
    except (IndexError, ValueError) as exc:
        raise ValueError(
            f"{filename}: could not read PID, x_fin, y_fin, and z_fin."
        ) from exc

    ids_position = arr[:, 0].astype(np.int64)
    x_fin_values = arr[:, 1].astype(float)
    y_fin_values = arr[:, 2].astype(float)
    z_fin_values = arr[:, 3].astype(float)

    return (
        ids_position,
        x_fin_values,
        y_fin_values,
        z_fin_values,
        tuple(coordinate_columns),
    )


def calculate_theta_fin_degrees(x_fin, y_fin, z_fin, mode):
    """Calculate the requested final-position angle in degrees."""
    cylindrical_radius = np.hypot(x_fin, y_fin)

    if mode == "axis":
        theta = np.arctan2(cylindrical_radius, np.abs(z_fin))
    elif mode == "polar":
        theta = np.arctan2(cylindrical_radius, z_fin)
    elif mode == "latitude":
        theta = np.arctan2(z_fin, cylindrical_radius)
    else:
        raise ValueError(f"Unknown theta mode: {mode}")

    theta_degrees = np.degrees(theta)
    invalid_position = (
        ~np.isfinite(x_fin)
        | ~np.isfinite(y_fin)
        | ~np.isfinite(z_fin)
        | ((cylindrical_radius == 0.0) & (z_fin == 0.0))
    )
    theta_degrees[invalid_position] = np.nan
    return theta_degrees


fn_stat = "stat_traj.dat"
ids_stat, mass_stat_raw, ye_fin_raw, s_fin_raw, vr_fin_raw, ye_5GK_raw, s_5GK_raw, bernoulli_raw, r_fin_raw,_ = stat.read_stat_traj(fn_stat)

ids_time_fin_stat, time_fin_raw, time_fin_column = read_time_fin_from_stat(
    fn_stat,
    column_override=args.time_fin_column,
)

ids_time_rfl_stat, time_rfl_raw, time_rfl_column = read_time_rfl_from_stat(
    fn_stat,
    column_override=args.time_rfl_column,
)

(
    ids_position_stat,
    x_fin_raw,
    y_fin_raw,
    z_fin_raw,
    position_columns,
) = read_final_positions_from_stat(fn_stat)

fn_dye = "weak_freezeout_dye.dat"

ids = np.loadtxt(fn_dye,comments='#', usecols = 0, dtype=int)
data = np.loadtxt(fn_dye,comments='#', usecols = (1,2,9,10,11,12,13,14,15,16), dtype = float)

mass = data[:,0]
t_FO_dye = data[:,1]
rho_FO_dye = data[:,2]
T_FO_dye = data[:,3]/MeV_to_K
Ye_FO_dye = data[:,4]
Rec_FO_dye = data[:,5]
Rpc_FO_dye = data[:,6]
t_exp_FO_dye = data[:,7]
t_weak_FO_dye = data[:,8]
dye_after_FO_dye = data[:,9]

fn_time = "weak_freezeout_timescale.dat"

ids_time = np.loadtxt(fn_time, comments='#', usecols=0, dtype=int)
data = np.loadtxt(fn_time,comments='#', usecols = (1,2,9,10,11,12,13,14,15,16), dtype = float)

if not np.array_equal(ids_time, ids):
    raise ValueError(
        "The particle-ID order differs between weak_freezeout_dye.dat "
        "and weak_freezeout_timescale.dat. Align the timescale data by PID "
        "before using ids to identify histogram entries."
    )

t_FO_time = data[:,1]
rho_FO_time = data[:,2]
T_FO_time = data[:,3]/MeV_to_K
Ye_FO_time = data[:,4]
Rec_FO_time = data[:,5]
Rpc_FO_time = data[:,6]
t_exp_FO_time = data[:,7]
t_weak_FO_time = data[:,8]
dye_after_FO_time = data[:,9]


fn_dng = "weak_freezeout_dng.dat"

ids_dng = np.loadtxt(fn_dng, comments='#', usecols=0, dtype=int)
data = np.loadtxt(fn_dng,comments='#', usecols = (1,2,9,10,11,12,13,14,15,16), dtype = float)

if not np.array_equal(ids_dng, ids):
    raise ValueError(
        "The particle-ID order differs between weak_freezeout_dye.dat "
        "and weak_freezeout_dng.dat. Align the dng data by PID "
        "before using ids to identify histogram entries."
    )

t_FO_dng = data[:,1]
rho_FO_dng = data[:,2]
T_FO_dng = data[:,3]/MeV_to_K
Ye_FO_dng = data[:,4]
Rec_FO_dng = data[:,5]
Rpc_FO_dng = data[:,6]
t_exp_FO_dng = data[:,7]
t_weak_FO_dng = data[:,8]
dye_after_FO_dng = data[:,9]



fn_dnv = "weak_freezeout_dnv.dat"

ids_dnv = np.loadtxt(fn_dnv, comments='#', usecols=0, dtype=int)
data = np.loadtxt(fn_dnv,comments='#', usecols = (1,2,9,10,11,12,13,14,15,16), dtype = float)

if not np.array_equal(ids_dnv, ids):
    raise ValueError(
        "The particle-ID order differs between weak_freezeout_dye.dat "
        "and weak_freezeout_dnv.dat. Align the dnv data by PID "
        "before using ids to identify histogram entries."
    )

t_FO_dnv = data[:,1]
rho_FO_dnv = data[:,2]
T_FO_dnv = data[:,3]/MeV_to_K
Ye_FO_dnv = data[:,4]
Rec_FO_dnv = data[:,5]
Rpc_FO_dnv = data[:,6]
t_exp_FO_dnv = data[:,7]
t_weak_FO_dnv = data[:,8]
dye_after_FO_dnv = data[:,9]



mass_stat, ye_fin, s_fin, vr_fin, ye_5GK, s_5GK, bernoulli, r_fin = stat.align_to_reference_ids(
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
    label=fn_stat,
)

(time_fin,) = stat.align_to_reference_ids(
    ids,
    ids_time_fin_stat,
    time_fin_raw,
    label=f"{fn_stat}:time_fin",
)

(time_rfl,) = stat.align_to_reference_ids(
    ids,
    ids_time_rfl_stat,
    time_rfl_raw,
    label=f"{fn_stat}:time_rfl",
)

x_fin, y_fin, z_fin = stat.align_to_reference_ids(
    ids,
    ids_position_stat,
    x_fin_raw,
    y_fin_raw,
    z_fin_raw,
    label=f"{fn_stat}:final_position",
)

theta_fin = calculate_theta_fin_degrees(
    x_fin,
    y_fin,
    z_fin,
    mode=args.theta_mode,
)

weights = mass

# Keep the existing Bernoulli-based ejecta definition.
minimum_mask = (
    np.isfinite(weights)
    & (weights > 0.0)
    & np.isfinite(bernoulli)
    & (bernoulli < 0.0)
)

ejecta_mask = minimum_mask.copy()

ejecta_mask = stat.add_range_condition(
    ejecta_mask,
    time_fin,
    args.time_fin_min,
    args.time_fin_max,
    open_upper=args.open_upper,
)

ejecta_mask = stat.add_range_condition(
    ejecta_mask,
    time_rfl,
    args.time_rfl_min,
    args.time_rfl_max,
    open_upper=args.open_upper,
)

ejecta_mask = stat.add_range_condition(
    ejecta_mask,
    theta_fin,
    args.theta_fin_min,
    args.theta_fin_max,
    open_upper=args.open_upper,
)

mass_ej_tot = np.sum(mass[minimum_mask])
if not np.isfinite(mass_ej_tot) or mass_ej_tot <= 0.0:
    raise ValueError(
        "No positive tracer mass remains after the Bernoulli ejecta cut."
    )


def range_selection_text(quantity, lower, upper, unit=""):
    """Return a compact description of one active numerical range cut."""
    if lower is None and upper is None:
        return ""

    unit_suffix = f" {unit}" if unit else ""
    if lower is not None and upper is not None:
        upper_relation = "<" if args.open_upper else "<="
        return (
            f"{lower:.10g} <= {quantity} {upper_relation} "
            f"{upper:.10g}{unit_suffix}"
        )
    if lower is not None:
        return f"{quantity} >= {lower:.10g}{unit_suffix}"

    upper_relation = "<" if args.open_upper else "<="
    return f"{quantity} {upper_relation} {upper:.10g}{unit_suffix}"


time_fin_cut_text = range_selection_text(
    "time_fin",
    args.time_fin_min,
    args.time_fin_max,
    unit="s",
)
time_rfl_cut_text = range_selection_text(
    "time_rfl",
    args.time_rfl_min,
    args.time_rfl_max,
    unit="s",
)
theta_cut_text = range_selection_text(
    f"theta_fin[{args.theta_mode}]",
    args.theta_fin_min,
    args.theta_fin_max,
    unit="deg",
)

active_primary_cuts = [
    text
    for text in (time_fin_cut_text, time_rfl_cut_text, theta_cut_text)
    if text
]
primary_cut_suffix = (
    " with " + ", ".join(active_primary_cuts)
    if active_primary_cuts
    else ""
)



def _latex_range(name, lower, upper, unit=""):
    """Return one LaTeX-formatted active cut, or None if inactive."""
    if lower is None and upper is None:
        return None

    upper_relation = "<" if args.open_upper else r"\leq"
    if lower is not None and upper is not None:
        return rf"${lower:.3g}\leq {name}{upper_relation}{upper:.3g}{unit}$"
    if lower is not None:
        return rf"${name}\geq {lower:.3g}{unit}$"
    return rf"${name}{upper_relation}{upper:.3g}{unit}$"


def build_cut_annotation():
    """Build a multiline description of all cuts used for common_mask."""
    lines = [r"$\mathrm{Cuts:}$", r"$h_{\rm Bernoulli}<0$"]

    candidate_lines = (
        _latex_range(
            r"t_{\rm fin}",
            args.time_fin_min,
            args.time_fin_max,
            r"\,\mathrm{s}",
        ),
        _latex_range(
            r"t_{\rm rfl}",
            args.time_rfl_min,
            args.time_rfl_max,
            r"\,\mathrm{s}",
        ),
        _latex_range(
            r"Y_{\rm e}(5\,\mathrm{GK})",
            args.ye_min,
            args.ye_max,
        ),
        _latex_range(
            r"s(5\,\mathrm{GK})/k_{\rm B}",
            args.s_min,
            args.s_max,
        ),
        _latex_range(
            r"v^r_{\rm fin}",
            args.vr_min,
            args.vr_max,
            r"\,\mathrm{cm\,s^{-1}}",
        ),
        _latex_range(
            r"r_{\rm fin}",
            args.r_fin_min,
            args.r_fin_max,
            r"\,\mathrm{cm}",
        ),
        _latex_range(
            {
                "axis": r"\theta_{\rm axis}",
                "polar": r"\theta_{\rm polar}",
                "latitude": r"\theta_{\rm lat}",
            }[args.theta_mode],
            args.theta_fin_min,
            args.theta_fin_max,
            r"^{\circ}",
        ),
    )

    lines.extend(line for line in candidate_lines if line is not None)
    return "\n".join(lines)


def add_cut_annotation(ax, fontsize=11, loc="upper left"):
    """Draw the active-cut description inside an axes."""
    locations = {
        "upper left": (0.02, 0.98, "left", "top"),
        "upper right": (0.98, 0.98, "right", "top"),
        "lower left": (0.02, 0.02, "left", "bottom"),
        "lower right": (0.98, 0.02, "right", "bottom"),
    }
    if loc not in locations:
        raise ValueError(
            f"Unknown cut-label location {loc!r}; choose from {tuple(locations)}."
        )
    x, y, ha, va = locations[loc]

    ax.text(
        x,
        y,
        build_cut_annotation(),
        transform=ax.transAxes,
        ha=ha,
        va=va,
        fontsize=fontsize,
        zorder=20,
        bbox=dict(
            boxstyle="round,pad=0.35",
            facecolor="white",
            edgecolor="0.7",
            alpha=0.85,
        ),
    )

selected_ejecta_id_file = (
    f"{args.output_prefix}_selected_ids.txt"
    if active_primary_cuts
    else f"{args.output_prefix}_ejecta_ids.txt"
)
np.savetxt(
    selected_ejecta_id_file,
    ids[ejecta_mask],
    fmt="%d",
    header=(
        "PID selected by bernoulli < 0"
        + primary_cut_suffix
    ),
)
print(f"Wrote selected ejecta IDs: {selected_ejecta_id_file}")

common_mask = ejecta_mask.copy()
common_mask = stat.add_range_condition(
    common_mask, ye_5GK, args.ye_min, args.ye_max,
    open_upper=args.open_upper
)
common_mask = stat.add_range_condition(
    common_mask, s_5GK, args.s_min, args.s_max,
    open_upper=args.open_upper
)

common_mask = stat.add_range_condition(
    common_mask, vr_fin, args.vr_min, args.vr_max,
    open_upper=args.open_upper
)

common_mask = stat.add_range_condition(
    common_mask, r_fin, args.r_fin_min, args.r_fin_max,
    open_upper=args.open_upper
)

common_selected_id_file = (
    f"{args.output_prefix}_common_selected_ids.txt"
)

np.savetxt(
    common_selected_id_file,
    ids[common_mask],
    fmt="%d",
    header=(
        "PID selected by all active cuts used for the T_FO histograms"
    ),
)

print(
    f"Wrote fully selected IDs: {common_selected_id_file}"
)


# Use the freeze-out time belonging to each definition.
mask_dye = stat.add_range_condition(
    common_mask, t_FO_dye, None, None,
    open_upper=args.open_upper
)
mask_time = stat.add_range_condition(
    common_mask, t_FO_time, None, None,
    open_upper=args.open_upper
)
mask_dng = stat.add_range_condition(
    common_mask, t_FO_dng, None, None,
    open_upper=args.open_upper
)
mask_dnv = stat.add_range_condition(
    common_mask, t_FO_dnv, None, None,
    open_upper=args.open_upper
)


selection_description = "Bernoulli-selected ejecta" + primary_cut_suffix

print(
    f"{selection_description}: "
    f"n={np.count_nonzero(ejecta_mask)}, "
    f"M_ej={mass_ej_tot:.10e} g"
)
print(
    f"time_fin was read from 1-based column {time_fin_column + 1} "
    f"of {fn_stat}."
)
print(
    f"time_rfl was read from 1-based column {time_rfl_column + 1} "
    f"of {fn_stat}."
)
print(
    "Final position was read from 1-based columns "
    f"{position_columns[0] + 1}, {position_columns[1] + 1}, "
    f"and {position_columns[2] + 1} of {fn_stat}."
)

print(
    f"time_fin range: {args.time_fin_min}, {args.time_fin_max}"
)
print(
    f"time_rfl range: {args.time_rfl_min}, {args.time_rfl_max}"
)
print(
    f"theta_fin mode/range [deg]: {args.theta_mode}, "
    f"{args.theta_fin_min}, {args.theta_fin_max}"
)
print(
    f"Histogram selections: "
    f"dye n={np.count_nonzero(mask_dye)}, "
    f"time n={np.count_nonzero(mask_time)}"
    f"dng n={np.count_nonzero(mask_dng)}"
    f"dnv n={np.count_nonzero(mask_dnv)}"
)

# str1 = ""
# for pid in ids[ (base_mask) & (T_FO_time < 0.3) ]:
#     if len(str1)>0:
#         str1 += ","
#         pass
#     str1 += "traj_%08d.dat" % (pid)
#     pass
# print(str1)

bins_ye = np.linspace(0.005,0.605,61)

fig_yehist=plt.figure(figsize=(10, 10*0.625))
ax_yehist=fig_yehist.add_subplot(111)

# col = "k"
# label = r"$\Delta Y_\mathrm{e} (t>t_\mathrm{FO}) < 0.01$"
# ax_yehist.hist(Ye_FO_dye[mask_dye], weights=mass[mask_dye]/mass_ej_tot, bins=bins_ye, histtype="step", log=True, alpha=0.7, align='mid', lw=2.0, color=col, label=label)

# col = "r"
# label = r"$t_\mathrm{weak}>t_\mathrm{exp}$"
# ax_yehist.hist(Ye_FO_time[mask_time], weights=mass[mask_time]/mass_ej_tot, bins=bins_ye, histtype="step", log=True, alpha=0.7, align='mid', lw=2.0, color=col, label=label)

# col = "g"
# label = r"$\Delta N_\mathrm{gross} < 0.01$"
# ax_yehist.hist(Ye_FO_dng[mask_dng], weights=mass[mask_dng]/mass_ej_tot, bins=bins_ye, histtype="step", log=True, alpha=0.7, align='mid', lw=2.0, color=col, label=label)

# col = "y"
# label = r"$\Delta N_\mathrm{var} < 0.01$"
# ax_yehist.hist(Ye_FO_dnv[mask_dnv], weights=mass[mask_dnv]/mass_ej_tot, bins=bins_ye, histtype="step", log=True, alpha=0.7, align='mid', lw=2.0, color=col, label=label)

col = "0.3"
label = r"all $Y_\mathrm{e}(5\,\mathrm{GK})$"
ax_yehist.hist(
    ye_5GK,
    weights=mass / mass_ej_tot,
    bins=bins_ye,
    histtype='step',
    log=True,
    alpha=0.3,
    align='mid',
    lw=2.0,
    color=col,
    label=label,
)

col = "r"
label = r"Selected $Y_\mathrm{e}(5\,\mathrm{GK})$"
ax_yehist.hist(
    ye_5GK[mask_dye],
    weights=mass[mask_dye] / mass_ej_tot,
    bins=bins_ye,
    histtype='step',
    log=True,
    alpha=0.7,
    align='mid',
    lw=2.0,
    color=col,
    label=label,
)

fig_yehist.subplots_adjust(left=0.15,bottom=0.15,right=0.95, top=0.95)
ax_yehist.legend(fontsize=16,frameon=False,ncol=1)
ax_yehist.set_xlim(0.0, 0.60)
ax_yehist.set_ylim(1e-5, 1)
ax_yehist.set_yscale('log')
ax_yehist.tick_params(axis='both',width=1.5,length=6,which='major',top=True,right=True,pad=10.0)
ax_yehist.tick_params(axis='both',width=1.5,length=3,which='minor',top=True,right=True)
ax_yehist.set_xlabel(r"$Y_\mathrm{e}$")
ax_yehist.set_ylabel("$\\Delta M/M_{\\mathrm{ej}}$")

ax_yehist.xaxis.set_major_locator(ticker.MultipleLocator(0.1))
ax_yehist.xaxis.set_minor_locator(ticker.MultipleLocator(0.02))



bins_s = np.logspace(0.0,3.0,101)

fig_shist=plt.figure(figsize=(10, 10*0.625))
ax_shist=fig_shist.add_subplot(111)

# col = "k"
# label = r"$\Delta Y_\mathrm{e} (t>t_\mathrm{FO}) < 0.01$"
# ax_shist.hist(Ye_FO_dye[mask_dye], weights=mass[mask_dye]/mass_ej_tot, bins=bins_ye, histtype="step", log=True, alpha=0.7, align='mid', lw=2.0, color=col, label=label)

# col = "r"
# label = r"$t_\mathrm{weak}>t_\mathrm{exp}$"
# ax_shist.hist(Ye_FO_time[mask_time], weights=mass[mask_time]/mass_ej_tot, bins=bins_ye, histtype="step", log=True, alpha=0.7, align='mid', lw=2.0, color=col, label=label)

col = "0.3"
label = r"all $s(5\,\mathrm{GK})$"
ax_shist.hist(
    s_5GK,
    weights=mass / mass_ej_tot,
    bins=bins_s,
    histtype='step',
    log=True,
    alpha=0.3,
    align='mid',
    lw=2.0,
    color=col,
    label=label,
)

col = "b"
label = r"Selected $s(5\,\mathrm{GK})$"
ax_shist.hist(
    s_5GK[mask_dye],
    weights=mass[mask_dye] / mass_ej_tot,
    bins=bins_s,
    histtype='step',
    log=True,
    alpha=0.7,
    align='mid',
    lw=2.0,
    color=col,
    label=label,
)

fig_shist.subplots_adjust(left=0.15,bottom=0.15,right=0.95, top=0.95)
ax_shist.legend(fontsize=16,frameon=False,ncol=1)
ax_shist.set_xlim(3.0, 300.0)
ax_shist.set_ylim(1e-5, 1)
ax_shist.set_xscale('log')
ax_shist.set_yscale('log')
ax_shist.tick_params(axis='both',width=1.5,length=6,which='major',top=True,right=True,pad=10.0)
ax_shist.tick_params(axis='both',width=1.5,length=3,which='minor',top=True,right=True)
ax_shist.set_xlabel(r"$s/k$")
ax_shist.set_ylabel("$\\Delta M/M_{\\mathrm{ej}}$")

#ax_shist.xaxis.set_major_locator(ticker.MultipleLocator(0.05))
#ax_shist.xaxis.set_minor_locator(ticker.MultipleLocator(0.01))
# The figure is saved after the optional T_FO-selected components
# are added below.


bins_T = np.logspace(-1,1,51)

T_hist_dye, _ = np.histogram(
    T_FO_dye[mask_dye],
    bins=bins_T,
    weights=mass[mask_dye] / mass_ej_tot,
)

T_hist_time, _ = np.histogram(
    T_FO_time[mask_time],
    bins=bins_T,
    weights=mass[mask_time] / mass_ej_tot,
)

T_hist_dng, _ = np.histogram(
    T_FO_dng[mask_dng],
    bins=bins_T,
    weights=mass[mask_dng] / mass_ej_tot,
)

T_hist_dnv, _ = np.histogram(
    T_FO_dnv[mask_dnv],
    bins=bins_T,
    weights=mass[mask_dnv] / mass_ej_tot,
)

T_center = np.sqrt(bins_T[:-1] * bins_T[1:])

histogram_output = np.column_stack(
    (
        bins_T[:-1],
        bins_T[1:],
        T_center,
        T_hist_dye,
        T_hist_time,
        T_hist_dng,
        T_hist_dnv,
    )
)

np.savetxt(
    f"{args.output_prefix}_TFO_histograms.txt",
    histogram_output,
    header=(
        "T_lo[MeV] T_hi[MeV] T_center[MeV] "
        "dM_over_Mej_dye dM_over_Mej_time "
        "dM_over_Mej_dng dM_over_Mej_dnv"
    ),
)


def get_temperature_bin(target_temperature):
    """Return the histogram-bin index and edges containing target_temperature."""
    ibin = np.searchsorted(bins_T, target_temperature, side="right") - 1
    if ibin < 0 or ibin >= len(bins_T) - 1:
        raise ValueError(
            f"Inspection temperature {target_temperature:g} MeV is outside "
            f"[{bins_T[0]:g}, {bins_T[-1]:g}] MeV."
        )
    return ibin, bins_T[ibin], bins_T[ibin + 1]


def temperature_tag(value):
    """Return a compact filename-safe temperature tag."""
    return (
        f"{value:.8g}"
        .replace("-", "m")
        .replace("+", "")
        .replace(".", "p")
    )


def write_particles_in_temperature_selection(
    source_name,
    selection_name,
    selection_tag,
    T_lo,
    T_hi,
    bin_index,
    requested_temperature,
    selection_mask,
    temperature,
    freezeout_time,
    freezeout_rho,
    freezeout_ye,
    freezeout_rec,
    freezeout_rpc,
    freezeout_t_exp,
    freezeout_t_weak,
    dye_after_freezeout,
):
    """Write tracer information for one temperature interval."""
    if not np.isfinite(T_lo) or not np.isfinite(T_hi):
        raise ValueError(
            f"Temperature limits must be finite: {T_lo}, {T_hi}"
        )
    if T_lo >= T_hi:
        raise ValueError(
            "Temperature lower limit must be smaller than upper limit: "
            f"{T_lo:g} >= {T_hi:g} MeV."
        )

    in_range = stat.add_range_condition(
        selection_mask,
        temperature,
        T_lo,
        T_hi,
        open_upper=args.open_upper,
    )
    in_range &= np.isfinite(temperature)

    selected_ids = ids[in_range]
    selected_mass = mass[in_range]

    output_data = np.column_stack(
        (
            selected_ids,
            freezeout_time[in_range],
            freezeout_rho[in_range],
            temperature[in_range],
            freezeout_ye[in_range],
            freezeout_rec[in_range],
            freezeout_rpc[in_range],
            freezeout_t_exp[in_range],
            freezeout_t_weak[in_range],
            dye_after_freezeout[in_range],
            time_fin[in_range],
            time_rfl[in_range],
            x_fin[in_range],
            y_fin[in_range],
            z_fin[in_range],
            theta_fin[in_range],
            r_fin[in_range],
            vr_fin[in_range],
            ye_fin[in_range],
            s_fin[in_range],
            ye_5GK[in_range],
            s_5GK[in_range],
            bernoulli[in_range],
            selected_mass,
            selected_mass / mass_ej_tot,
        )
    )

    output_file = (
        f"{args.output_prefix}_{selection_tag}_particles.txt"
    )
    trajectory_file = (
        f"{args.output_prefix}_{selection_tag}_trajectories.txt"
    )

    upper_symbol_text = "<" if args.open_upper else "<="
    interval_text = (
        f"{T_lo:.10g} <= T_FO {upper_symbol_text} {T_hi:.10g} MeV"
    )

    if requested_temperature is None:
        request_text = "Requested exact range"
    else:
        request_text = (
            f"Requested temperature: {requested_temperature:.10g} MeV"
        )

    center = (
        np.sqrt(T_lo * T_hi)
        if T_lo > 0.0 and T_hi > 0.0
        else 0.5 * (T_lo + T_hi)
    )

    header = (
        f"Histogram source: {source_name}\n"
        f"Selection type: {selection_name}\n"
        f"{request_text}\n"
        f"Bin index: {bin_index}\n"
        f"Selected interval: {interval_text}\n"
        f"Representative temperature: {center:.10g} MeV\n"
        f"Number of particles: {selected_ids.size}\n"
        "Columns: "
        "pid  t_FO[s]  rho_FO[g/cm^3]  T_FO[MeV]  Ye_FO  "
        "Rec_FO[1/s]  Rpc_FO[1/s]  t_exp_FO[s]  t_weak_FO[s]  "
        "dYe_after_FO  time_fin[s]  time_rfl[s]  "
        "x_fin[cm]  y_fin[cm]  z_fin[cm]  "
        f"theta_fin_{args.theta_mode}[deg]  "
        "r_fin[cm]  vr_fin[cm/s]  "
        "Ye_fin  s_fin[kB/nuc]  Ye_5GK  s_5GK[kB/nuc]  "
        "bernoulli_column  mass[g]  mass/M_ej"
    )

    np.savetxt(
        output_file,
        output_data,
        fmt=["%d"] + ["%.10e"] * 24,
        header=header,
    )

    with open(trajectory_file, "w") as f:
        for pid in selected_ids:
            f.write(f"traj_{pid:08d}.dat\n")

    print(
        f"[{source_name}] {interval_text}: "
        f"n={selected_ids.size}, "
        f"mass/M_ej={np.sum(selected_mass) / mass_ej_tot:.8e}"
    )
    print(f"  tracer table: {output_file}")
    print(f"  trajectory list: {trajectory_file}")

    return {
        "source": source_name,
        "selection_name": selection_name,
        "target": requested_temperature,
        "ibin": bin_index,
        "T_lo": T_lo,
        "T_hi": T_hi,
        "T_center": center,
        "mask": in_range,
        "n": selected_ids.size,
        "mass": np.sum(selected_mass),
    }



def write_all_cut_selected_particles(
    source_name,
    selection_mask,
    temperature,
    freezeout_time,
    freezeout_rho,
    freezeout_ye,
    freezeout_rec,
    freezeout_rpc,
    freezeout_t_exp,
    freezeout_t_weak,
    dye_after_freezeout,
):
    """
    Write one particles table without applying a T_FO-bin/range cut.

    The output contains all tracers that pass selection_mask and have finite
    t_FO and T_FO for the requested freeze-out definition.
    """
    selected = np.asarray(selection_mask, dtype=bool).copy()
    selected &= np.isfinite(freezeout_time)
    selected &= np.isfinite(temperature)

    selected_ids = ids[selected]
    selected_mass = mass[selected]

    output_data = np.column_stack(
        (
            selected_ids,
            freezeout_time[selected],
            freezeout_rho[selected],
            temperature[selected],
            freezeout_ye[selected],
            freezeout_rec[selected],
            freezeout_rpc[selected],
            freezeout_t_exp[selected],
            freezeout_t_weak[selected],
            dye_after_freezeout[selected],
            time_fin[selected],
            time_rfl[selected],
            x_fin[selected],
            y_fin[selected],
            z_fin[selected],
            theta_fin[selected],
            r_fin[selected],
            vr_fin[selected],
            ye_fin[selected],
            s_fin[selected],
            ye_5GK[selected],
            s_5GK[selected],
            bernoulli[selected],
            selected_mass,
            selected_mass / mass_ej_tot,
        )
    )

    output_file = f"{args.output_prefix}_selected_{source_name}_particles.txt"
    trajectory_file = f"{args.output_prefix}_selected_{source_name}_trajectories.txt"

    header = (
        f"Freeze-out source: {source_name}\n"
        "Selection type: all active non-temperature cuts; no T_FO cut\n"
        f"Number of particles: {selected_ids.size}\n"
        "Columns: "
        "pid  t_FO[s]  rho_FO[g/cm^3]  T_FO[MeV]  Ye_FO  "
        "Rec_FO[1/s]  Rpc_FO[1/s]  t_exp_FO[s]  t_weak_FO[s]  "
        "dYe_after_FO  time_fin[s]  time_rfl[s]  "
        "x_fin[cm]  y_fin[cm]  z_fin[cm]  "
        f"theta_fin_{args.theta_mode}[deg]  "
        "r_fin[cm]  vr_fin[cm/s]  "
        "Ye_fin  s_fin[kB/nuc]  Ye_5GK  s_5GK[kB/nuc]  "
        "bernoulli_column  mass[g]  mass/M_ej"
    )

    np.savetxt(
        output_file,
        output_data,
        fmt=["%d"] + ["%.10e"] * 24,
        header=header,
    )

    with open(trajectory_file, "w") as stream:
        for pid in selected_ids:
            stream.write(f"traj_{pid:08d}.dat\n")

    print(
        f"[{source_name}] all active cuts, no T_FO cut: "
        f"n={selected_ids.size}, "
        f"mass/M_ej={np.sum(selected_mass) / mass_ej_tot:.8e}"
    )
    print(f"  tracer table: {output_file}")
    print(f"  trajectory list: {trajectory_file}")

    return {
        "source": source_name,
        "mask": selected,
        "n": selected_ids.size,
        "mass": np.sum(selected_mass),
        "output_file": output_file,
        "trajectory_file": trajectory_file,
    }

def weighted_summary(values, selected_weights):
    """Return the 15.865, 50, and 84.135 mass-weighted percentiles."""
    valid = (
        np.isfinite(values)
        & np.isfinite(selected_weights)
        & (selected_weights > 0.0)
    )
    if not np.any(valid):
        return np.array([np.nan, np.nan, np.nan])
    return stat.weighted_quantile(
        values[valid],
        selected_weights[valid],
        np.array([0.15865, 0.5, 0.84135]),
    )


def write_particles_in_temperature_bin(
    source_name,
    target_temperature,
    selection_mask,
    temperature,
    freezeout_time,
    freezeout_rho,
    freezeout_ye,
    freezeout_rec,
    freezeout_rpc,
    freezeout_t_exp,
    freezeout_t_weak,
    dye_after_freezeout,
):
    """Write tracer information for the histogram bin containing T."""
    ibin, T_lo, T_hi = get_temperature_bin(target_temperature)
    center = np.sqrt(T_lo * T_hi)

    return write_particles_in_temperature_selection(
        source_name=source_name,
        selection_name="histogram bin",
        selection_tag=(
            f"Tbin_{source_name}_{temperature_tag(center)}MeV"
        ),
        T_lo=T_lo,
        T_hi=T_hi,
        bin_index=ibin,
        requested_temperature=target_temperature,
        selection_mask=selection_mask,
        temperature=temperature,
        freezeout_time=freezeout_time,
        freezeout_rho=freezeout_rho,
        freezeout_ye=freezeout_ye,
        freezeout_rec=freezeout_rec,
        freezeout_rpc=freezeout_rpc,
        freezeout_t_exp=freezeout_t_exp,
        freezeout_t_weak=freezeout_t_weak,
        dye_after_freezeout=dye_after_freezeout,
    )


def write_particles_in_temperature_range(
    source_name,
    T_lo,
    T_hi,
    selection_mask,
    temperature,
    freezeout_time,
    freezeout_rho,
    freezeout_ye,
    freezeout_rec,
    freezeout_rpc,
    freezeout_t_exp,
    freezeout_t_weak,
    dye_after_freezeout,
):
    """Write tracer information for an exact user-specified range."""
    range_tag = (
        f"Trange_{source_name}_{temperature_tag(T_lo)}to"
        f"{temperature_tag(T_hi)}MeV"
    )

    return write_particles_in_temperature_selection(
        source_name=source_name,
        selection_name="exact user-specified range",
        selection_tag=range_tag,
        T_lo=T_lo,
        T_hi=T_hi,
        bin_index=-1,
        requested_temperature=None,
        selection_mask=selection_mask,
        temperature=temperature,
        freezeout_time=freezeout_time,
        freezeout_rho=freezeout_rho,
        freezeout_ye=freezeout_ye,
        freezeout_rec=freezeout_rec,
        freezeout_rpc=freezeout_rpc,
        freezeout_t_exp=freezeout_t_exp,
        freezeout_t_weak=freezeout_t_weak,
        dye_after_freezeout=dye_after_freezeout,
    )


def write_requested_temperature_groups(
    source_name,
    selection_mask,
    temperature,
    freezeout_time,
    freezeout_rho,
    freezeout_ye,
    freezeout_rec,
    freezeout_rpc,
    freezeout_t_exp,
    freezeout_t_weak,
    dye_after_freezeout,
):
    """Dispatch to exact-range mode or backward-compatible bin mode."""
    groups = []

    if args.inspect_T_ranges:
        for T_lo, T_hi in args.inspect_T_ranges:
            groups.append(
                write_particles_in_temperature_range(
                    source_name,
                    T_lo,
                    T_hi,
                    selection_mask,
                    temperature,
                    freezeout_time,
                    freezeout_rho,
                    freezeout_ye,
                    freezeout_rec,
                    freezeout_rpc,
                    freezeout_t_exp,
                    freezeout_t_weak,
                    dye_after_freezeout,
                )
            )
    else:
        for target_temperature in args.inspect_T:
            groups.append(
                write_particles_in_temperature_bin(
                    source_name,
                    target_temperature,
                    selection_mask,
                    temperature,
                    freezeout_time,
                    freezeout_rho,
                    freezeout_ye,
                    freezeout_rec,
                    freezeout_rpc,
                    freezeout_t_exp,
                    freezeout_t_weak,
                    dye_after_freezeout,
                )
            )

    return groups


def write_group_comparison(groups, source_arrays):
    """Write compact mass-weighted comparisons of selected T_FO groups."""
    if not groups:
        return

    overview_file = f"{args.output_prefix}_peak_groups_overview.txt"
    overview_rows = []
    for igroup, group in enumerate(groups, start=1):
        overview_rows.append(
            [
                igroup,
                group["ibin"],
                group["T_lo"],
                group["T_hi"],
                group["T_center"],
                group["n"],
                group["mass"],
                group["mass"] / mass_ej_tot,
            ]
        )
    np.savetxt(
        overview_file,
        np.asarray(overview_rows),
        fmt=["%d", "%d", "%.10e", "%.10e", "%.10e", "%d", "%.10e", "%.10e"],
        header=(
            "group  bin_index(-1=exact_range)  T_lo[MeV]  "
            "T_hi[MeV]  T_center[MeV]  n_tracer  mass[g]  mass/M_ej"
        ),
    )

    quantities_file = f"{args.output_prefix}_peak_groups_quantiles.txt"
    quantity_rows = []
    for igroup, group in enumerate(groups, start=1):
        group_mask = group["mask"]
        group_weights = mass[group_mask]
        for quantity_name, values in source_arrays.items():
            p16, p50, p84 = weighted_summary(values[group_mask], group_weights)
            quantity_rows.append(
                (
                    igroup,
                    group["T_center"],
                    quantity_name,
                    p16,
                    p50,
                    p84,
                )
            )

    with open(quantities_file, "w") as f:
        f.write(
            "# group  T_center[MeV]  quantity  "
            "p15p865_mass_weighted  p50_mass_weighted  "
            "p84p135_mass_weighted\n"
        )
        for row in quantity_rows:
            f.write(
                f"{row[0]:d} {row[1]:.10e} {row[2]:s} "
                f"{row[3]:.10e} {row[4]:.10e} {row[5]:.10e}\n"
            )

    print(f"Wrote peak-group overview: {overview_file}")
    print(f"Wrote peak-group quantiles: {quantities_file}")



# When no T_FO inspection is requested, still write a source-specific table
# containing every tracer that passes the active cuts.
if not args.inspect_T and not args.inspect_T_ranges:
    if args.inspect_source in ("dye", "both"):
        write_all_cut_selected_particles(
            "dye",
            mask_dye,
            T_FO_dye,
            t_FO_dye,
            rho_FO_dye,
            Ye_FO_dye,
            Rec_FO_dye,
            Rpc_FO_dye,
            t_exp_FO_dye,
            t_weak_FO_dye,
            dye_after_FO_dye,
        )

    if args.inspect_source in ("time", "both"):
        write_all_cut_selected_particles(
            "time",
            mask_time,
            T_FO_time,
            t_FO_time,
            rho_FO_time,
            Ye_FO_time,
            Rec_FO_time,
            Rpc_FO_time,
            t_exp_FO_time,
            t_weak_FO_time,
            dye_after_FO_time,
        )

groups_by_source = {}

if args.inspect_source in ("dye", "both"):
    dye_groups = write_requested_temperature_groups(
        "dye",
        mask_dye,
        T_FO_dye,
        t_FO_dye,
        rho_FO_dye,
        Ye_FO_dye,
        Rec_FO_dye,
        Rpc_FO_dye,
        t_exp_FO_dye,
        t_weak_FO_dye,
        dye_after_FO_dye,
    )
    if dye_groups:
        groups_by_source["dye"] = dye_groups

    tweak_over_texp_dye = np.full_like(t_weak_FO_dye, np.nan)
    np.divide(
        t_weak_FO_dye,
        t_exp_FO_dye,
        out=tweak_over_texp_dye,
        where=np.isfinite(t_exp_FO_dye) & (t_exp_FO_dye != 0.0),
    )

    write_group_comparison(
        dye_groups,
        {
            "t_FO_s": t_FO_dye,
            "rho_FO_gcm3": rho_FO_dye,
            "T_FO_MeV": T_FO_dye,
            "Ye_FO": Ye_FO_dye,
            "Rec_FO_1ps": Rec_FO_dye,
            "Rpc_FO_1ps": Rpc_FO_dye,
            "t_exp_FO_s": t_exp_FO_dye,
            "t_weak_FO_s": t_weak_FO_dye,
            "tweak_over_texp": tweak_over_texp_dye,
            "dYe_after_FO": dye_after_FO_dye,
            "time_fin_s": time_fin,
            "time_rfl_s": time_rfl,
            "x_fin_cm": x_fin,
            "y_fin_cm": y_fin,
            "z_fin_cm": z_fin,
            f"theta_fin_{args.theta_mode}_deg": theta_fin,
            "r_fin_cm": r_fin,
            "vr_fin_cmps": vr_fin,
            "Ye_fin": ye_fin,
            "s_fin_kBpnuc": s_fin,
            "Ye_5GK": ye_5GK,
            "s_5GK_kBpnuc": s_5GK,
            "bernoulli_column": bernoulli,
        },
    )

if args.inspect_source in ("time", "both"):
    time_groups = write_requested_temperature_groups(
        "time",
        mask_time,
        T_FO_time,
        t_FO_time,
        rho_FO_time,
        Ye_FO_time,
        Rec_FO_time,
        Rpc_FO_time,
        t_exp_FO_time,
        t_weak_FO_time,
        dye_after_FO_time,
    )
    if time_groups:
        groups_by_source["time"] = time_groups


# Use the same colors for
#   * T_FO-selected components in the Ye_FO histogram, and
#   * highlight_groups in make_figure_selected_5GK().
#
# Each mask is in the same PID order as ids because ids_ref=ids is supplied
# to make_figure_selected_5GK() below.
highlight_groups = []
upper_relation = "<" if args.open_upper else r"\leq"
group_colors = plt.rcParams["axes.prop_cycle"].by_key().get(
    "color",
    [f"C{i}" for i in range(10)],
)
icolor = 0

freezeout_ye_by_source = {
    "dye": Ye_FO_dye,
    "time": Ye_FO_time,
}

for source_name, source_groups in groups_by_source.items():
    for group in source_groups:
        color = group_colors[icolor % len(group_colors)]
        icolor += 1

        label = (
            rf"{source_name}: "
            rf"${group['T_lo']:.2f}\leq T_{{\rm FO}}"
            rf"{upper_relation}{group['T_hi']:.2f}\ {{\rm MeV}}$"
        )

        highlight_groups.append(
            {
                "mask": group["mask"],
                "label": label,
                "color": color,
            }
        )

        selected = group["mask"]

        # For exact user-specified T_FO ranges, compare the corresponding
        # Ye(5 GK) distribution. For the older histogram-bin mode, keep
        # the previous behavior and show Ye at freeze-out.
        if args.inspect_T_ranges:
            selected_ye_for_hist = ye_5GK[selected]
        else:
            selected_ye_for_hist = freezeout_ye_by_source[source_name][selected]

        ax_yehist.hist(
            selected_ye_for_hist,
            weights=mass[selected] / mass_ej_tot,
            bins=bins_ye,
            histtype="step",
            log=True,
            alpha=0.95,
            align="mid",
            lw=2.5,
            color=color,
            label=label,
        )

# Rebuild the legend after adding the optional T_FO-selected components.
ax_yehist.legend(
    fontsize=12,
    frameon=False,
    ncol=1,
    loc="best",
)
add_cut_annotation(
    ax_yehist,
    fontsize=args.cut_label_fontsize,
    loc=args.cut_label_loc,
)
fn_fig = f"Ye_FO_{args.output_prefix}.pdf"
fig_yehist.savefig(fn_fig)


ax_shist.legend(
    fontsize=12,
    frameon=False,
    ncol=1,
    loc="best",
)
add_cut_annotation(
    ax_shist,
    fontsize=args.cut_label_fontsize,
    loc=args.cut_label_loc,
)
fn_fig = f"s_FO_{args.output_prefix}.pdf"
fig_shist.savefig(fn_fig)

# Older alternative:
# stat.make_figure_selected(
#     ye_5GK,
#     s_5GK,
#     vr_fin,
#     args.ye_min,
#     args.ye_max,
#     args.s_min,
#     args.s_max,
#     args.vr_min,
#     args.vr_max,
#     args.output_prefix,
#     highlight_groups=highlight_groups,
# )
make_figure_selected_5GK(
    stat_file=fn_stat,
    prefix=args.output_prefix,
    ids_ref=ids,
    background_mask=common_mask,
    highlight_groups=highlight_groups,
    s_range=(3.0, 300.0),
    texp_range=(1e-4, 1e1),
    annotation_text=build_cut_annotation(),
    annotation_fontsize=args.cut_label_fontsize,
    annotation_loc=args.cut_label_loc,
)

fig_TFOhist_dye=plt.figure(figsize=(10, 10*0.625))
ax_TFOhist_dye=fig_TFOhist_dye.add_subplot(111)

col = "k"
label = r"$\Delta Y_\mathrm{e} (t>t_\mathrm{FO}) < 0.01$"
ax_TFOhist_dye.hist(T_FO_dye[minimum_mask], weights=mass[minimum_mask]/mass_ej_tot, bins=bins_T, histtype="step", log=True, alpha=0.4, align='mid', ls="dashed", lw=2.0, color=col, label="")
ax_TFOhist_dye.hist(T_FO_dye[mask_dye], weights=mass[mask_dye]/mass_ej_tot, bins=bins_T, histtype="step", log=True, alpha=0.7, align='mid', lw=1.0, color=col, label=label)

fig_TFOhist_dye.subplots_adjust(left=0.15,bottom=0.15,right=0.95, top=0.95)
ax_TFOhist_dye.legend(fontsize=16,frameon=False,ncol=1)
ax_TFOhist_dye.tick_params(axis='both',width=1.5,length=6,which='major',top=True,right=True,pad=10.0)
ax_TFOhist_dye.tick_params(axis='both',width=1.5,length=3,which='minor',top=True,right=True)

ax_TFOhist_dye.set_xlim(0.1, 10.0)
ax_TFOhist_dye.set_xscale('log')
ax_TFOhist_dye.set_xlabel(r"$T$ (MeV)")

ax_TFOhist_dye.set_ylim(1e-5, 1)
ax_TFOhist_dye.set_yscale('log')
ax_TFOhist_dye.set_ylabel("$\\Delta M/M_{\\mathrm{ej}}$")

#ax_TFOhist_dye.xaxis.set_major_locator(ticker.MultipleLocator(0.05))
#ax_TFOhist_dye.xaxis.set_minor_locator(ticker.MultipleLocator(0.01))
add_cut_annotation(
    ax_TFOhist_dye,
    fontsize=args.cut_label_fontsize,
    loc=args.cut_label_loc,
)
fn_fig = f"T_FO_dye_{args.output_prefix}.pdf"
fig_TFOhist_dye.savefig(fn_fig)

fig_TFOhist_time=plt.figure(figsize=(10, 10*0.625))
ax_TFOhist_time=fig_TFOhist_time.add_subplot(111)

col = "r"
label = r"$t_\mathrm{weak}>t_\mathrm{exp}$"
ax_TFOhist_time.hist(T_FO_time[minimum_mask], weights=mass[minimum_mask]/mass_ej_tot, bins=bins_T, histtype="step", log=True, alpha=0.4, align='mid', ls="dashed", lw=2.0, color=col, label="")
ax_TFOhist_time.hist(T_FO_time[mask_time], weights=mass[mask_time]/mass_ej_tot, bins=bins_T, histtype="step", log=True, alpha=0.7, align='mid', lw=1.0, color=col, label=label)

fig_TFOhist_time.subplots_adjust(left=0.15,bottom=0.15,right=0.95, top=0.95)
ax_TFOhist_time.legend(fontsize=16,frameon=False,ncol=1)
ax_TFOhist_time.tick_params(axis='both',width=1.5,length=6,which='major',top=True,right=True,pad=10.0)
ax_TFOhist_time.tick_params(axis='both',width=1.5,length=3,which='minor',top=True,right=True)

ax_TFOhist_time.set_xlim(0.1, 10.0)
ax_TFOhist_time.set_xscale('log')
ax_TFOhist_time.set_xlabel(r"$T$ (MeV)")

ax_TFOhist_time.set_ylim(1e-5, 1)
ax_TFOhist_time.set_yscale('log')
ax_TFOhist_time.set_ylabel("$\\Delta M/M_{\\mathrm{ej}}$")

#ax_TFOhist_time.xaxis.set_major_locator(ticker.MultipleLocator(0.05))
#ax_TFOhist_time.xaxis.set_minor_locator(ticker.MultipleLocator(0.01))
add_cut_annotation(
    ax_TFOhist_time,
    fontsize=args.cut_label_fontsize,
    loc=args.cut_label_loc,
)
fn_fig = f"T_FO_time_{args.output_prefix}.pdf"
fig_TFOhist_time.savefig(fn_fig)

fig_TFOhist_dng=plt.figure(figsize=(10, 10*0.625))
ax_TFOhist_dng=fig_TFOhist_dng.add_subplot(111)

col = "g"
label = r"$\Delta N_\mathrm{gross} < 0.01$"
ax_TFOhist_dng.hist(T_FO_dng[minimum_mask], weights=mass[minimum_mask]/mass_ej_tot, bins=bins_T, histtype="step", log=True, alpha=0.4, align='mid', ls="dashed", lw=2.0, color=col, label="")
ax_TFOhist_dng.hist(T_FO_dng[mask_dng], weights=mass[mask_dng]/mass_ej_tot, bins=bins_T, histtype="step", log=True, alpha=0.7, align='mid', lw=1.0, color=col, label=label)

fig_TFOhist_dng.subplots_adjust(left=0.15,bottom=0.15,right=0.95, top=0.95)
ax_TFOhist_dng.legend(fontsize=16,frameon=False,ncol=1)
ax_TFOhist_dng.tick_params(axis='both',width=1.5,length=6,which='major',top=True,right=True,pad=10.0)
ax_TFOhist_dng.tick_params(axis='both',width=1.5,length=3,which='minor',top=True,right=True)

ax_TFOhist_dng.set_xlim(0.1, 10.0)
ax_TFOhist_dng.set_xscale('log')
ax_TFOhist_dng.set_xlabel(r"$T$ (MeV)")

ax_TFOhist_dng.set_ylim(1e-5, 1)
ax_TFOhist_dng.set_yscale('log')
ax_TFOhist_dng.set_ylabel("$\\Delta M/M_{\\mathrm{ej}}$")

#ax_TFOhist_dng.xaxis.set_major_locator(ticker.MultipleLocator(0.05))
#ax_TFOhist_dng.xaxis.set_minor_locator(ticker.MultipleLocator(0.01))
add_cut_annotation(
    ax_TFOhist_dng,
    fontsize=args.cut_label_fontsize,
    loc=args.cut_label_loc,
)
fn_fig = f"T_FO_dng_{args.output_prefix}.pdf"
fig_TFOhist_dng.savefig(fn_fig)

fig_TFOhist_dnv=plt.figure(figsize=(10, 10*0.625))
ax_TFOhist_dnv=fig_TFOhist_dnv.add_subplot(111)

col = "b"
label = r"$\Delta N_\mathrm{var} < 0.01$"
ax_TFOhist_dnv.hist(T_FO_dnv[minimum_mask], weights=mass[minimum_mask]/mass_ej_tot, bins=bins_T, histtype="step", log=True, alpha=0.4, align='mid', ls="dashed", lw=2.0, color=col, label="")
ax_TFOhist_dnv.hist(T_FO_dnv[mask_dnv], weights=mass[mask_dnv]/mass_ej_tot, bins=bins_T, histtype="step", log=True, alpha=0.7, align='mid', lw=1.0, color=col, label=label)

fig_TFOhist_dnv.subplots_adjust(left=0.15,bottom=0.15,right=0.95, top=0.95)
ax_TFOhist_dnv.legend(fontsize=16,frameon=False,ncol=1)
ax_TFOhist_dnv.tick_params(axis='both',width=1.5,length=6,which='major',top=True,right=True,pad=10.0)
ax_TFOhist_dnv.tick_params(axis='both',width=1.5,length=3,which='minor',top=True,right=True)

ax_TFOhist_dnv.set_xlim(0.1, 10.0)
ax_TFOhist_dnv.set_xscale('log')
ax_TFOhist_dnv.set_xlabel(r"$T$ (MeV)")

ax_TFOhist_dnv.set_ylim(1e-5, 1)
ax_TFOhist_dnv.set_yscale('log')
ax_TFOhist_dnv.set_ylabel("$\\Delta M/M_{\\mathrm{ej}}$")

#ax_TFOhist_dnv.xaxis.set_major_locator(ticker.MultipleLocator(0.05))
#ax_TFOhist_dnv.xaxis.set_minor_locator(ticker.MultipleLocator(0.01))
add_cut_annotation(
    ax_TFOhist_dnv,
    fontsize=args.cut_label_fontsize,
    loc=args.cut_label_loc,
)
fn_fig = f"T_FO_dnv_{args.output_prefix}.pdf"
fig_TFOhist_dnv.savefig(fn_fig)


plt.show()
