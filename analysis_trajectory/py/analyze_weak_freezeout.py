#!/usr/bin/env python3
"""Select ejecta and compare four weak-freeze-out definitions."""

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


MEV_TO_K = 1.16e10

# stat_traj.dat columns (0-based). The analysis format has 36 columns.
COL_ID = 0
COL_MASS = 1
COL_TIME_FIN = 2
COL_X_FIN = 3
COL_Y_FIN = 4
COL_Z_FIN = 5
COL_VR_FIN = 6
COL_S_FIN = 10
COL_YE_FIN = 11
COL_S_5GK = 13
COL_YE_5GK = 14
COL_TEXP_5GK = 15
COL_BERNOULLI = 18
COL_TIME_RFL = 35
N_STAT_COLUMNS = 36


def read_freezeout_file(filename):
    """Read the columns shared by the four weak_freezeout_*.dat files."""
    data = np.loadtxt(
        filename,
        comments="#",
        usecols=(0, 1, 2, 9, 10, 11, 12, 13, 14, 15, 16),
        ndmin=2,
    )
    return data[:, 0].astype(np.int64), {
        "mass": data[:, 1],
        "time": data[:, 2],
        "rho": data[:, 3],
        "temperature": data[:, 4] / MEV_TO_K,
        "ye": data[:, 5],
        "rec": data[:, 6],
        "rpc": data[:, 7],
        "texp": data[:, 8],
        "tweak": data[:, 9],
        "dye_after": data[:, 10],
    }


parser = argparse.ArgumentParser(
    description="Select weak-freeze-out tracers from the fixed 36-column format."
)
parser.add_argument("--time-rfl-min", "--time_rfl_min", type=float, default=0.0)
parser.add_argument("--time-rfl-max", "--time_rfl_max", type=float)
parser.add_argument("--theta-min", "--theta-fin-min", "--theta_fin_min", type=float)
parser.add_argument("--theta-max", "--theta-fin-max", "--theta_fin_max", type=float)
parser.add_argument("--ye-5gk-min", type=float)
parser.add_argument("--ye-5gk-max", type=float)
parser.add_argument("--s-5gk-min", type=float)
parser.add_argument("--s-5gk-max", type=float)
parser.add_argument("--vr-min", type=float)
parser.add_argument("--vr-max", type=float)
parser.add_argument("--r-fin-min", type=float)
parser.add_argument("--r-fin-max", type=float)
parser.add_argument("--time-fin-min", "--time_fin_min", type=float)
parser.add_argument("--time-fin-max", "--time_fin_max", type=float)
parser.add_argument(
    "--open-upper",
    action="store_true",
    help="Use < rather than <= for every upper bound.",
)
parser.add_argument(
    "--inspect-T",
    type=float,
    nargs="+",
    default=[],
    help="Select the histogram bin containing each temperature [MeV].",
)
parser.add_argument(
    "--inspect-T-range",
    dest="inspect_T_ranges",
    type=float,
    nargs=2,
    action="append",
    default=[],
    metavar=("T_MIN", "T_MAX"),
)
parser.add_argument(
    "--inspect-source",
    choices=("dye", "time", "both"),
    default="dye",
)
parser.add_argument("-o", "--output-prefix", default="weakfo")
args = parser.parse_args()

if args.inspect_T and args.inspect_T_ranges:
    parser.error("--inspect-T and --inspect-T-range cannot be used together")

range_arguments = (
    ("time-rfl", args.time_rfl_min, args.time_rfl_max),
    ("theta", args.theta_min, args.theta_max),
    ("ye-5gk", args.ye_5gk_min, args.ye_5gk_max),
    ("s-5gk", args.s_5gk_min, args.s_5gk_max),
    ("vr", args.vr_min, args.vr_max),
    ("r-fin", args.r_fin_min, args.r_fin_max),
    ("time-fin", args.time_fin_min, args.time_fin_max),
)
for name, lower, upper in range_arguments:
    if lower is not None and not np.isfinite(lower):
        parser.error(f"--{name}-min must be finite")
    if upper is not None and not np.isfinite(upper):
        parser.error(f"--{name}-max must be finite")
    if lower is not None and upper is not None and lower >= upper:
        parser.error(f"--{name}-min must be smaller than --{name}-max")

for name, value in (("theta-min", args.theta_min), ("theta-max", args.theta_max)):
    if value is not None and not 0.0 <= value <= 90.0:
        parser.error(f"--{name} must be between 0 and 90 degrees")

prefix = Path(args.output_prefix)
if prefix.suffix.lower() in (".txt", ".dat"):
    prefix = prefix.with_suffix("")
prefix.parent.mkdir(parents=True, exist_ok=True)
stem = prefix.name


# Read the fixed stat_traj.dat format explicitly.
stat = np.loadtxt("stat_traj.dat", comments="#", ndmin=2)
if stat.shape[1] < N_STAT_COLUMNS:
    raise ValueError(
        f"stat_traj.dat: expected {N_STAT_COLUMNS} columns, got {stat.shape[1]}"
    )

ids_stat = stat[:, COL_ID].astype(np.int64)
mass_stat = stat[:, COL_MASS]
time_fin_stat = stat[:, COL_TIME_FIN]
x_fin_stat = stat[:, COL_X_FIN]
y_fin_stat = stat[:, COL_Y_FIN]
z_fin_stat = stat[:, COL_Z_FIN]
vr_fin_stat = stat[:, COL_VR_FIN]
s_fin_stat = stat[:, COL_S_FIN]
ye_fin_stat = stat[:, COL_YE_FIN]
s_5gk_stat = stat[:, COL_S_5GK]
ye_5gk_stat = stat[:, COL_YE_5GK]
texp_5gk_stat = stat[:, COL_TEXP_5GK]
bernoulli_stat = stat[:, COL_BERNOULLI]
time_rfl_stat = stat[:, COL_TIME_RFL]


# Read the four freeze-out definitions. Their PID order is expected to match.
ids, freezeout_dye = read_freezeout_file("weak_freezeout_dye.dat")
ids_time, freezeout_time = read_freezeout_file("weak_freezeout_timescale.dat")
ids_dng, freezeout_dng = read_freezeout_file("weak_freezeout_dng.dat")
ids_dnv, freezeout_dnv = read_freezeout_file("weak_freezeout_dnv.dat")

for filename, ids_other in (
    ("weak_freezeout_timescale.dat", ids_time),
    ("weak_freezeout_dng.dat", ids_dng),
    ("weak_freezeout_dnv.dat", ids_dnv),
):
    if not np.array_equal(ids_other, ids):
        raise ValueError(f"PID order in {filename} differs from weak_freezeout_dye.dat")


# Align stat_traj.dat to the PID order of the freeze-out tables.
stat_index = {int(pid): index for index, pid in enumerate(ids_stat)}
if len(stat_index) != len(ids_stat):
    raise ValueError("stat_traj.dat contains duplicate PIDs")
missing_ids = [int(pid) for pid in ids if int(pid) not in stat_index]
if missing_ids:
    raise ValueError(
        f"{len(missing_ids)} freeze-out PIDs are absent from stat_traj.dat; "
        f"first missing PID: {missing_ids[0]}"
    )
order = np.fromiter(
    (stat_index[int(pid)] for pid in ids), dtype=np.int64, count=len(ids)
)

mass_stat = mass_stat[order]
time_fin = time_fin_stat[order]
x_fin = x_fin_stat[order]
y_fin = y_fin_stat[order]
z_fin = z_fin_stat[order]
vr_fin = vr_fin_stat[order]
s_fin = s_fin_stat[order]
ye_fin = ye_fin_stat[order]
s_5gk = s_5gk_stat[order]
ye_5gk = ye_5gk_stat[order]
texp_5gk = texp_5gk_stat[order]
bernoulli = bernoulli_stat[order]
time_rfl = time_rfl_stat[order]

r_fin = np.sqrt(x_fin**2 + y_fin**2 + z_fin**2)
cylindrical_radius = np.hypot(x_fin, y_fin)
theta = np.degrees(np.arctan2(cylindrical_radius, np.abs(z_fin)))
invalid_position = (
    ~np.isfinite(x_fin)
    | ~np.isfinite(y_fin)
    | ~np.isfinite(z_fin)
    | ((cylindrical_radius == 0.0) & (z_fin == 0.0))
)
theta[invalid_position] = np.nan

# Keep the historical weights from weak_freezeout_dye.dat. They agree with
# stat_traj.dat for the production files and are the original T_FO weights.
mass = freezeout_dye["mass"]


# Bernoulli ejecta and the primary geometric/time cuts.
base_mask = (
    np.isfinite(mass)
    & (mass > 0.0)
    & np.isfinite(bernoulli)
    & (bernoulli < 0.0)
)
ejecta_mask = base_mask.copy()

if args.time_rfl_min is not None:
    ejecta_mask &= np.isfinite(time_rfl) & (time_rfl >= args.time_rfl_min)
if args.time_rfl_max is not None:
    if args.open_upper:
        ejecta_mask &= np.isfinite(time_rfl) & (time_rfl < args.time_rfl_max)
    else:
        ejecta_mask &= np.isfinite(time_rfl) & (time_rfl <= args.time_rfl_max)

if args.theta_min is not None:
    ejecta_mask &= np.isfinite(theta) & (theta >= args.theta_min)
if args.theta_max is not None:
    if args.open_upper:
        ejecta_mask &= np.isfinite(theta) & (theta < args.theta_max)
    else:
        ejecta_mask &= np.isfinite(theta) & (theta <= args.theta_max)

if args.time_fin_min is not None:
    ejecta_mask &= np.isfinite(time_fin) & (time_fin >= args.time_fin_min)
if args.time_fin_max is not None:
    if args.open_upper:
        ejecta_mask &= np.isfinite(time_fin) & (time_fin < args.time_fin_max)
    else:
        ejecta_mask &= np.isfinite(time_fin) & (time_fin <= args.time_fin_max)


# The remaining cuts are applied to the common population used by all four
# freeze-out definitions.
common_mask = ejecta_mask.copy()

if args.ye_5gk_min is not None:
    common_mask &= np.isfinite(ye_5gk) & (ye_5gk >= args.ye_5gk_min)
if args.ye_5gk_max is not None:
    if args.open_upper:
        common_mask &= np.isfinite(ye_5gk) & (ye_5gk < args.ye_5gk_max)
    else:
        common_mask &= np.isfinite(ye_5gk) & (ye_5gk <= args.ye_5gk_max)

if args.s_5gk_min is not None:
    common_mask &= np.isfinite(s_5gk) & (s_5gk >= args.s_5gk_min)
if args.s_5gk_max is not None:
    if args.open_upper:
        common_mask &= np.isfinite(s_5gk) & (s_5gk < args.s_5gk_max)
    else:
        common_mask &= np.isfinite(s_5gk) & (s_5gk <= args.s_5gk_max)

if args.vr_min is not None:
    common_mask &= np.isfinite(vr_fin) & (vr_fin >= args.vr_min)
if args.vr_max is not None:
    if args.open_upper:
        common_mask &= np.isfinite(vr_fin) & (vr_fin < args.vr_max)
    else:
        common_mask &= np.isfinite(vr_fin) & (vr_fin <= args.vr_max)

if args.r_fin_min is not None:
    common_mask &= np.isfinite(r_fin) & (r_fin >= args.r_fin_min)
if args.r_fin_max is not None:
    if args.open_upper:
        common_mask &= np.isfinite(r_fin) & (r_fin < args.r_fin_max)
    else:
        common_mask &= np.isfinite(r_fin) & (r_fin <= args.r_fin_max)


mass_ejecta = np.sum(mass[base_mask])
mass_selected = np.sum(mass[common_mask])
if not np.isfinite(mass_ejecta) or mass_ejecta <= 0.0:
    raise ValueError("No positive mass remains after the Bernoulli ejecta cut")

np.savetxt(
    prefix.parent / f"{stem}_selected_ids.txt",
    ids[ejecta_mask],
    fmt="%d",
    header="PID selected by Bernoulli criterion and primary time/angle cuts",
)
np.savetxt(
    prefix.parent / f"{stem}_common_selected_ids.txt",
    ids[common_mask],
    fmt="%d",
    header="PID selected by all active cuts used for the T_FO histograms",
)


# Mass-weighted freeze-out temperature histograms. The normalization stays
# the total positive-mass Bernoulli ejecta mass, as in the original script.
freezeout = {
    "dye": freezeout_dye,
    "time": freezeout_time,
    "dng": freezeout_dng,
    "dnv": freezeout_dnv,
}
definition_masks = {
    name: common_mask & np.isfinite(values["time"])
    for name, values in freezeout.items()
}

temperature_edges = np.logspace(-1.0, 1.0, 51)
temperature_center = np.sqrt(temperature_edges[:-1] * temperature_edges[1:])
temperature_histograms = {}
for name, values in freezeout.items():
    selected = definition_masks[name]
    temperature_histograms[name], _ = np.histogram(
        values["temperature"][selected],
        bins=temperature_edges,
        weights=mass[selected] / mass_ejecta,
    )

np.savetxt(
    prefix.parent / f"{stem}_TFO_histograms.txt",
    np.column_stack(
        (
            temperature_edges[:-1],
            temperature_edges[1:],
            temperature_center,
            temperature_histograms["dye"],
            temperature_histograms["time"],
            temperature_histograms["dng"],
            temperature_histograms["dnv"],
        )
    ),
    header=(
        "T_lo[MeV] T_hi[MeV] T_center[MeV] "
        "dM_over_Mej_dye dM_over_Mej_time "
        "dM_over_Mej_dng dM_over_Mej_dnv"
    ),
)


# --inspect-T and --inspect-T-range optionally narrow the particle tables.
def inspected_temperature_mask(temperatures):
    if not args.inspect_T and not args.inspect_T_ranges:
        return np.ones(len(temperatures), dtype=bool)

    selected = np.zeros(len(temperatures), dtype=bool)
    for target in args.inspect_T:
        index = np.searchsorted(temperature_edges, target, side="right") - 1
        if index < 0 or index >= len(temperature_edges) - 1:
            raise ValueError(
                f"Inspection temperature {target:g} MeV is outside "
                f"[{temperature_edges[0]:g}, {temperature_edges[-1]:g}] MeV"
            )
        selected |= (
            np.isfinite(temperatures)
            & (temperatures >= temperature_edges[index])
            & (temperatures < temperature_edges[index + 1])
        )

    for lower, upper in args.inspect_T_ranges:
        if not np.isfinite(lower) or not np.isfinite(upper) or lower >= upper:
            raise ValueError("Each --inspect-T-range must have finite T_MIN < T_MAX")
        in_range = np.isfinite(temperatures) & (temperatures >= lower)
        if args.open_upper:
            in_range &= temperatures < upper
        else:
            in_range &= temperatures <= upper
        selected |= in_range

    return selected


def write_selected_particles(source, selected):
    """Write the downstream PID table and trajectory filename list."""
    values = freezeout[source]
    table = np.column_stack(
        (
            ids[selected],
            values["time"][selected],
            values["rho"][selected],
            values["temperature"][selected],
            values["ye"][selected],
            values["rec"][selected],
            values["rpc"][selected],
            values["texp"][selected],
            values["tweak"][selected],
            values["dye_after"][selected],
            time_fin[selected],
            time_rfl[selected],
            x_fin[selected],
            y_fin[selected],
            z_fin[selected],
            theta[selected],
            r_fin[selected],
            vr_fin[selected],
            ye_fin[selected],
            s_fin[selected],
            ye_5gk[selected],
            s_5gk[selected],
            bernoulli[selected],
            mass[selected],
            mass[selected] / mass_ejecta,
        )
    )
    if args.inspect_T:
        selection_text = "active cuts and requested T_FO histogram bins"
    elif args.inspect_T_ranges:
        selection_text = "active cuts and requested T_FO ranges"
    else:
        selection_text = "all active non-temperature cuts; no T_FO cut"

    particle_file = prefix.parent / f"{stem}_selected_{source}_particles.txt"
    trajectory_file = prefix.parent / f"{stem}_selected_{source}_trajectories.txt"
    header = (
        f"Freeze-out source: {source}\n"
        f"Selection type: {selection_text}\n"
        f"Number of particles: {np.count_nonzero(selected)}\n"
        "Columns: pid  t_FO[s]  rho_FO[g/cm^3]  T_FO[MeV]  Ye_FO  "
        "Rec_FO[1/s]  Rpc_FO[1/s]  t_exp_FO[s]  t_weak_FO[s]  "
        "dYe_after_FO  time_fin[s]  time_rfl[s]  x_fin[cm]  y_fin[cm]  "
        "z_fin[cm]  theta_fin_axis[deg]  r_fin[cm]  vr_fin[cm/s]  "
        "Ye_fin  s_fin[kB/nuc]  Ye_5GK  s_5GK[kB/nuc]  "
        "bernoulli_column  mass[g]  mass/M_ej"
    )
    np.savetxt(
        particle_file,
        table,
        fmt=["%d"] + ["%.10e"] * 24,
        header=header,
    )
    np.savetxt(
        trajectory_file,
        np.asarray([f"traj_{pid:08d}.dat" for pid in ids[selected]]),
        fmt="%s",
    )
    print(
        f"[{source}] selected particles: n={np.count_nonzero(selected)}, "
        f"mass={np.sum(mass[selected]):.10e} g"
    )
    print(f"  particle table: {particle_file}")
    print(f"  trajectory list: {trajectory_file}")


sources_to_write = (
    ("dye", "time") if args.inspect_source == "both" else (args.inspect_source,)
)
for source in sources_to_write:
    selected = definition_masks[source] & inspected_temperature_mask(
        freezeout[source]["temperature"]
    )
    write_selected_particles(source, selected)


# Plotting is kept separate from the selection flow. Matplotlib mathtext is
# used so the script does not depend on an external LaTeX installation.
plt.rcParams.update(
    {
        "font.size": 14.0,
        "font.family": "DejaVu Serif",
        "axes.linewidth": 1.2,
        "xtick.direction": "in",
        "ytick.direction": "in",
        "text.usetex": False,
    }
)

cut_lines = ["Bernoulli < 0"]
for label, lower, upper in (
    ("time_rfl [s]", args.time_rfl_min, args.time_rfl_max),
    ("theta [deg]", args.theta_min, args.theta_max),
    ("time_fin [s]", args.time_fin_min, args.time_fin_max),
    ("Ye(5 GK)", args.ye_5gk_min, args.ye_5gk_max),
    ("s(5 GK)", args.s_5gk_min, args.s_5gk_max),
    ("vr_fin", args.vr_min, args.vr_max),
    ("r_fin", args.r_fin_min, args.r_fin_max),
):
    if lower is not None and upper is not None:
        relation = "<" if args.open_upper else "<="
        cut_lines.append(f"{lower:g} <= {label} {relation} {upper:g}")
    elif lower is not None:
        cut_lines.append(f"{label} >= {lower:g}")
    elif upper is not None:
        relation = "<" if args.open_upper else "<="
        cut_lines.append(f"{label} {relation} {upper:g}")
cut_text = "\n".join(cut_lines)


def add_cut_text(ax):
    ax.text(
        0.02,
        0.98,
        cut_text,
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=9,
        bbox={"facecolor": "white", "edgecolor": "0.7", "alpha": 0.85},
    )


definition_labels = {
    "dye": r"$\Delta Y_e < 0.01$",
    "time": r"$t_{\rm weak} > t_{\rm exp}$",
    "dng": r"$\Delta N_{\rm gross} < 0.01$",
    "dnv": r"$\Delta N_{\rm var} < 0.01$",
}
definition_colors = {
    "dye": "k",
    "time": "tab:red",
    "dng": "tab:green",
    "dnv": "tab:blue",
}

for name, values in freezeout.items():
    available = base_mask & np.isfinite(values["time"])
    selected = definition_masks[name]
    fig, ax = plt.subplots(figsize=(10.0, 6.25))
    ax.hist(
        values["temperature"][available],
        bins=temperature_edges,
        weights=mass[available] / mass_ejecta,
        histtype="step",
        linestyle="--",
        linewidth=1.5,
        color=definition_colors[name],
        alpha=0.5,
    )
    ax.hist(
        values["temperature"][selected],
        bins=temperature_edges,
        weights=mass[selected] / mass_ejecta,
        histtype="step",
        linewidth=2.0,
        color=definition_colors[name],
        label=definition_labels[name],
    )
    ax.set(xscale="log", yscale="log", xlim=(0.1, 10.0), ylim=(1.0e-5, 1.0))
    ax.set_xlabel(r"$T_{\rm FO}$ [MeV]")
    ax.set_ylabel(r"$\Delta M/M_{\rm ej}$")
    ax.legend(frameon=False)
    add_cut_text(ax)
    fig.tight_layout()
    fig.savefig(prefix.parent / f"T_FO_{name}_{stem}.pdf")
    plt.close(fig)


# Ye(5 GK) histogram with the common selected population overlaid.
valid_ye = np.isfinite(ye_5gk) & np.isfinite(mass) & (mass > 0.0)
fig, ax = plt.subplots(figsize=(10.0, 6.25))
ax.hist(
    ye_5gk[valid_ye],
    bins=np.linspace(0.005, 0.605, 61),
    weights=mass[valid_ye] / mass_ejecta,
    histtype="step",
    linewidth=1.5,
    color="0.4",
    alpha=0.5,
    label="all tracers",
)
ax.hist(
    ye_5gk[common_mask],
    bins=np.linspace(0.005, 0.605, 61),
    weights=mass[common_mask] / mass_ejecta,
    histtype="step",
    linewidth=2.0,
    color="tab:red",
    label="selected",
)
ax.set(xlim=(0.0, 0.6), yscale="log", ylim=(1.0e-5, 1.0))
ax.set_xlabel(r"$Y_e(5\,\mathrm{GK})$")
ax.set_ylabel(r"$\Delta M/M_{\rm ej}$")
ax.legend(frameon=False)
add_cut_text(ax)
fig.tight_layout()
fig.savefig(prefix.parent / f"Ye_FO_{stem}.pdf")
plt.close(fig)


def plot_scatter(x, y, xlabel, ylabel, filename, xlog=False, ylog=False):
    valid = np.isfinite(x) & np.isfinite(y)
    if xlog:
        valid &= x > 0.0
    if ylog:
        valid &= y > 0.0
    selected = valid & common_mask

    fig, ax = plt.subplots(figsize=(10.0, 6.25))
    ax.scatter(x[valid], y[valid], s=3, c="0.75", rasterized=True, label="all")
    ax.scatter(
        x[selected],
        y[selected],
        s=8,
        c="tab:red",
        edgecolors="none",
        rasterized=True,
        label="selected",
    )
    if xlog:
        ax.set_xscale("log")
    if ylog:
        ax.set_yscale("log")
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.legend(frameon=False, markerscale=2)
    add_cut_text(ax)
    fig.tight_layout()
    fig.savefig(prefix.parent / filename, dpi=180)
    plt.close(fig)


plot_scatter(
    s_5gk,
    ye_5gk,
    r"$s(5\,\mathrm{GK})/k_{\rm B}$",
    r"$Y_e(5\,\mathrm{GK})$",
    f"Ye-S_5GK-{stem}.png",
    xlog=True,
)
plot_scatter(
    texp_5gk,
    ye_5gk,
    r"$t_{\rm exp}(5\,\mathrm{GK})$ [s]",
    r"$Y_e(5\,\mathrm{GK})$",
    f"Ye-texp_5GK-{stem}.png",
    xlog=True,
)
plot_scatter(
    s_5gk,
    texp_5gk,
    r"$s(5\,\mathrm{GK})/k_{\rm B}$",
    r"$t_{\rm exp}(5\,\mathrm{GK})$ [s]",
    f"S-texp_5GK-{stem}.png",
    xlog=True,
    ylog=True,
)


print("time_rfl column: stat_traj.dat column 36 (0-based index 35, t_last_rfl)")
print(f"Bernoulli ejecta: n={np.count_nonzero(base_mask)}, mass={mass_ejecta:.10e} g")
print(
    f"Selected ejecta: n={np.count_nonzero(common_mask)}, "
    f"mass={mass_selected:.10e} g, mass/M_ej={mass_selected / mass_ejecta:.10e}"
)
print(
    "Histogram selections: "
    + ", ".join(
        f"{name} n={np.count_nonzero(mask)}"
        for name, mask in definition_masks.items()
    )
)
