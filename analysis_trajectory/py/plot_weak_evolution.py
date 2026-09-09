import numpy as np
import sys
import os

import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as colors
import matplotlib.ticker as ticker

import argparse

import tracer_analysis_utils as percentile
import re

def read_fortran_float(value):
    """
    Fortran出力に現れることがある
        2.364369-121
        1.234567+123
    を
        2.364369E-121
        1.234567E+123
    として読む。
    """
    if isinstance(value, bytes):
        value = value.decode("ascii")

    value = value.strip()

    try:
        return float(value)
    except ValueError:
        corrected = re.sub(
            r"^([+-]?(?:\d+(?:\.\d*)?|\.\d+))([+-]\d+)$",
            r"\1E\2",
            value,
        )

        if corrected == value:
            raise ValueError(f"数値に変換できません: {value!r}")

        return float(corrected)

MeV_to_GK = 11.6

show_in_MeV = True
show_in_GK  = False

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
    'xtick.direction'     : 'in'     ,
    'ytick.major.size'    : 6        ,
    'ytick.major.width'   : 1.5      ,
    'ytick.minor.size'    : 3        ,
    'ytick.minor.width'   : 1.5      ,
    'ytick.labelsize'     : 24.0     ,
    'ytick.direction'     : 'in'     ,
    'xtick.major.pad'     : 2        ,
    'xtick.minor.pad'     : 2        ,
    'ytick.major.pad'     : 2        ,
    'ytick.minor.pad'     : 2        ,
    'axes.linewidth'      : 1.5      ,
    'text.usetex'         : True    }

plt.rcParams.update(params)

prefix_default = "percentiles_selected_by_final_quantities"
parser = argparse.ArgumentParser()
parser.add_argument("--prefix", type=str, required=False, default=prefix_default, help="prefix of the files")
parser.add_argument('--plot_raw_traj', required=False, action='store_true')

parser.add_argument(
    "--tfo-hist-file",
    type=str,
    default=None,
    help=(
        "T_FO histogram table produced by "
        "plot_hist_weak_compare_dye_peaks_range_5GK_timefin_theta_annot.py"
    ),
)

args = parser.parse_args()
prefix = args.prefix

tfo_histograms = None
if args.tfo_hist_file is not None:
    tfo_data = np.loadtxt(
        args.tfo_hist_file,
        comments="#",
        ndmin=2,
    )

    T_lo = tfo_data[:, 0]
    T_hi = tfo_data[:, 1]
    T_edges = np.concatenate((T_lo, [T_hi[-1]]))

    def normalize_histogram(hist):
        hist = np.asarray(hist, dtype=float)
        valid = np.isfinite(hist) & (hist >= 0.0)
        result = np.zeros_like(hist)
        total = np.sum(hist[valid])

        if total > 0.0:
            result[valid] = hist[valid] / total

        return result

    tfo_histograms = {
        "edges": T_edges,
        "dye": normalize_histogram(tfo_data[:, 3]),
        "time": normalize_histogram(tfo_data[:, 4]),
        "dng": normalize_histogram(tfo_data[:, 5]),
        "dnv": normalize_histogram(tfo_data[:, 6]),
    }

def add_tfo_histograms(ax, histograms):
    if histograms is None:
        return None

    ax_tfo = ax.twinx()
    edges = histograms["edges"]

    ax_tfo.stairs(
        histograms["dye"],
        edges,
        color="k",
        lw=1.8,
        alpha=0.75,
        label=r"$\Delta Y_{\rm e}<0.01$",
    )
    ax_tfo.stairs(
        histograms["time"],
        edges,
        color="tab:red",
        lw=1.4,
        alpha=0.65,
        label=r"$t_{\rm weak}>t_{\rm exp}$",
    )
    ax_tfo.stairs(
        histograms["dng"],
        edges,
        color="tab:green",
        lw=1.4,
        alpha=0.65,
        label=r"$\Delta N_{\rm gross}<0.01$",
    )
    ax_tfo.stairs(
        histograms["dnv"],
        edges,
        color="tab:blue",
        lw=1.4,
        alpha=0.65,
        label=r"$\Delta N_{\rm var}<0.01$",
    )

    hist_max = max(
        np.max(histograms["dye"]),
        np.max(histograms["time"]),
        np.max(histograms["dng"]),
        np.max(histograms["dnv"]),
    )
    ax_tfo.set_ylim(
        0.0,
        1.15 * hist_max if hist_max > 0.0 else 1.0,
    )

    ax_tfo.set_yticks([])
    ax_tfo.set_ylabel("")
    ax_tfo.set_ylim(0.0,1.0)
    ax_tfo.tick_params(
        axis="y",
        which="both",
        right=False,
        labelright=False,
    )
    ax_tfo.spines["right"].set_visible(False)

    return ax_tfo

# fn = "median_sfho1215_trap_reduced.dat"

fn = prefix + "_ye.txt"
temp = np.loadtxt(fn,usecols=0)
ye_lw = np.loadtxt(fn,usecols=1)
ye_med = np.loadtxt(fn,usecols=2)
ye_rw = np.loadtxt(fn,usecols=3)

fn = prefix + "_yecap.txt"
yecap_lw = np.loadtxt(fn,usecols=1)
yecap_med = np.loadtxt(fn,usecols=2)
yecap_rw = np.loadtxt(fn,usecols=3)

fn = prefix + "_texp.txt"
texp_lw = np.loadtxt(fn,usecols=1)
texp_med = np.loadtxt(fn,usecols=2)
texp_rw = np.loadtxt(fn,usecols=3)

fn = prefix + "_caprate.txt"
rcap_lw = np.loadtxt(fn,usecols=1)
rcap_med = np.loadtxt(fn,usecols=2)
rcap_rw = np.loadtxt(fn,usecols=3)

fn = prefix + "_eta.txt"
eta_lw = np.loadtxt(fn,usecols=1)
eta_med = np.loadtxt(fn,usecols=2)
eta_rw = np.loadtxt(fn,usecols=3)

fn = prefix + "_r.txt"
r_lw = np.loadtxt(fn,usecols=1)
r_med = np.loadtxt(fn,usecols=2)
r_rw = np.loadtxt(fn,usecols=3)

fn = prefix + "_vr.txt"
vr_lw = np.loadtxt(fn,usecols=1)
vr_med = np.loadtxt(fn,usecols=2)
vr_rw = np.loadtxt(fn,usecols=3)

fn = prefix + "_time.txt"
time_lw = np.loadtxt(fn,usecols=1)
time_med = np.loadtxt(fn,usecols=2)
time_rw = np.loadtxt(fn,usecols=3)

fn = prefix + "_avtexp.txt"
avtexp_lw = np.loadtxt(fn,usecols=1)
avtexp_med = np.loadtxt(fn,usecols=2)
avtexp_rw = np.loadtxt(fn,usecols=3)

# rabs_lw = np.loadtxt(fn,usecols=10)
# rabs_med = np.loadtxt(fn,usecols=11)
# rabs_rw = np.loadtxt(fn,usecols=12)


fn = "traj_time"
temp_list, ids, _, time_temp = percentile.read_temperature_table(fn)
fn = "traj_Ye"
_, _, _, ye_temp = percentile.read_temperature_table(fn)
fn = "traj_yecap"
_, _, _, yecap_temp = percentile.read_temperature_table(fn)
fn = "traj_vel"
_, _, _, vr_temp = percentile.read_temperature_table(fn)



id_to_index = {int(pid): idx for idx, pid in enumerate(ids)}

# pid = 1701
# print("")
# idx = id_to_index[pid]
# print(idx)
# for it in range(len(time_temp[idx])):
#     print(time_temp[idx,it], temp_list[it]*11.6)
#     pass

# list_ids = [1, 10, 30, 100, 300, 1000, 3000, 10000, 30000]
# list_ids = ids[::4500]

fn = prefix + "_ids.txt"
ids_selected = np.loadtxt(fn, comments="#", dtype=int)
# print(ids_selected)
list_ids = ids_selected[::500]
#list_ids = [ids_selected[0]]
if not args.plot_raw_traj:
    list_ids = []
    pass

temp_GK = temp*MeV_to_GK

if show_in_MeV:
    temp_label = "$kT$ (MeV)"
    temp_show = temp
    temp_min = 0.3

    temp_max = 3.0
    temp_ticks_pos = [3, 1, 0.3]
    temp_ticks_label = ["3", "1", "0.3"]

    temp_max = 10.0
    temp_ticks_pos = [10, 3, 1, 0.3]
    temp_ticks_label = ["10","3", "1", "0.3"]

elif show_in_GK:
    temp_show = temp_GK
    temp_max = 30.0
    temp_min = 3.0
    temp_label = "$T$ ($10^9$K)"
    temp_ticks_pos = [30, 10, 3]
    temp_ticks_label = ["30", "10", "3"]
    pass

size_fig = 10.0
aratio = 0.6


fig_ye = plt.figure(figsize=(size_fig, size_fig*aratio))
ax_ye = fig_ye.add_subplot(111)

color = "k"
ax_ye.fill_between(temp_show, ye_lw, ye_rw, ls="solid", color=color, alpha=0.2, label="")
ax_ye.plot(temp_show, ye_med, ls="solid", marker="none", color=color, alpha=0.7, label="")

color = "r"
ax_ye.fill_between(temp_show, yecap_lw, yecap_rw, ls="solid", color=color, alpha=0.2, label="")
ax_ye.plot(temp_show, yecap_med, ls="solid", marker="none", color=color, alpha=0.7, label="$Y_\\mathrm{e}^\\mathrm{eq,cap}$")


ax_tfo_ye = add_tfo_histograms(
    ax_ye,
    tfo_histograms,
)


ax_ye.set_xlim(temp_max,temp_min)
ax_ye.set_xscale("log")
ax_ye.set_xlabel(temp_label)

ax_ye.set_ylim(0.0,0.55)
ax_ye.set_yscale("linear")
ax_ye.set_ylabel("$Y_\\mathrm{e}$")

fig_ye.subplots_adjust(left=0.15,bottom=0.15,right=0.95, top=0.95)
ax_ye.legend(fontsize=18,frameon=False,ncols=1,loc='lower right')
ax_ye.tick_params(axis='both', which='both', top=True, right=True, pad=10)

ax_ye.set_xticks(temp_ticks_pos)
ax_ye.set_xticklabels(temp_ticks_label)

ax_ye.xaxis.set_minor_formatter(ticker.NullFormatter())

ax_ye.yaxis.set_major_locator(ticker.MultipleLocator(0.1))
ax_ye.yaxis.set_minor_locator(ticker.MultipleLocator(0.02))


fig_ye.savefig("Ye-%s.pdf" % (prefix))


fig_yediff = plt.figure(figsize=(size_fig, size_fig*aratio))
ax_yediff = fig_yediff.add_subplot(111)

color = "k"
ax_yediff.plot(temp_show, ye_med-yecap_med, ls="solid", marker="none", color=color, alpha=0.7, label="")

ax_yediff.grid(lw=1.0, ls="dotted", color="grey", alpha=0.5)

ax_yediff.set_xlim(temp_max,temp_min)
ax_yediff.set_xscale("log")
ax_yediff.set_xlabel(temp_label)

ax_yediff.set_ylim(-0.3,0.3)
ax_yediff.set_yscale("linear")
ax_yediff.set_ylabel("$Y_\\mathrm{e}-Y_\\mathrm{e}^\\mathrm{eq,cap}$")

fig_yediff.subplots_adjust(left=0.15,bottom=0.15,right=0.95, top=0.95)
ax_yediff.legend(fontsize=18,frameon=False,ncols=1,loc='lower right')
ax_yediff.tick_params(axis='both', which='both', top=True, right=True, pad=10)

ax_yediff.set_xticks(temp_ticks_pos)
ax_yediff.set_xticklabels(temp_ticks_label)

ax_yediff.xaxis.set_minor_formatter(ticker.NullFormatter())

#ax_yediff.yaxis.set_major_locator(ticker.MultipleLocator(0.1))
#ax_yediff.yaxis.set_minor_locator(ticker.MultipleLocator(0.02))


fig_yediff.savefig("Yediff-%s.pdf" % (prefix))


fig_tscale = plt.figure(figsize=(size_fig, size_fig*aratio))
ax_tscale = fig_tscale.add_subplot(111)

color = "k"
ax_tscale.fill_between(temp_show, texp_lw, texp_rw, ls="solid", color=color, alpha=0.2, label="")
ax_tscale.plot(temp_show, texp_med, ls="solid", marker="none", color=color, alpha=0.7, label="Expansion")

# color = "g"
# ax_tscale.fill_between(temp_show, avtexp_lw, avtexp_rw, ls="solid", color=color, alpha=0.2, label="")
# ax_tscale.plot(temp_show, avtexp_med, ls="solid", marker="none", color=color, alpha=0.7, label="Expansion2")

color = "r"
ax_tscale.fill_between(temp_show, 1.0/rcap_lw, 1.0/rcap_rw, ls="solid", color=color, alpha=0.2, label="")
ax_tscale.plot(temp_show, 1.0/rcap_med, ls="solid", marker="none", color=color, alpha=0.7, label="$e^-e^+$ capture")


rcap_Tlimit =  np.log(2.0)/6146.0 * (1.0+3.0*1.26**2) * 24.0 * (temp/0.511)**5
rcap_elimit =  np.log(2.0)/6146.0 * (1.0+3.0*1.26**2) * eta_med**5/5.0* (temp/0.511)**5

#ax_tscale.plot(temp_show, 1.0/rcap_Tlimit, ls="dashed", marker="none", color=color, alpha=0.7, label="$\\eta_e\\ll 1$, $T/m_e\\gg1$ limit")
#ax_tscale.plot(temp_show, 1.0/rcap_elimit, ls="dotted", marker="none", color=color, alpha=0.7, label="$\\eta_e\\gg 1$, $T\\eta_e/m_e\\gg1$ limit")

# dlogr_med = np.diff(np.log(r_med))
# dlogr_lw = np.diff(np.log(r_lw))
# dlogr_rw = np.diff(np.log(r_rw))
# dt    = np.diff(time_med)
# dlogr_dt_med = dlogr_med/dt; dlogr_dt_med = np.append(dlogr_dt_med, dlogr_dt_med[-1])
# dlogr_dt_lw = dlogr_lw/dt; dlogr_dt_lw = np.append(dlogr_dt_lw, dlogr_dt_lw[-1])
# dlogr_dt_rw = dlogr_rw/dt; dlogr_dt_rw = np.append(dlogr_dt_rw, dlogr_dt_rw[-1])
# ax_tscale.plot(temp_show, 1.0/dlogr_dt_med, ls="solid", color="b", alpha=0.8, label="")
# ax_tscale.fill_between(temp_show, 1.0/dlogr_dt_lw, 1.0/dlogr_dt_rw, ls="solid", color="b", alpha=0.2, label="")

# color = "b"
# ax_tscale.fill_between(temp_show, 1.0/rabs_lw, 1.0/rabs_rw, ls="solid", color=color, alpha=0.2, label="")
# ax_tscale.plot(temp_show, 1.0/rabs_med, ls="solid", marker="none", color=color, alpha=0.7, label="$\\nu$ absorption")

ax_tfo_tscale = add_tfo_histograms(
    ax_tscale,
    tfo_histograms,
)

ax_tscale.set_xlim(temp_max,temp_min)
ax_tscale.set_xscale("log")
ax_tscale.set_xlabel(temp_label)

ax_tscale.set_ylim(1e-4, 1e1)
ax_tscale.set_yscale("log")
ax_tscale.set_ylabel("Timescale (s)")

fig_tscale.subplots_adjust(left=0.15,bottom=0.15,right=0.95, top=0.95)
ax_tscale.legend(fontsize=18,frameon=False,ncols=1,loc='upper left')
ax_tscale.tick_params(axis='both', which='both', top=True, right=True, pad=10)

ax_tscale.set_xticks(temp_ticks_pos)
ax_tscale.set_xticklabels(temp_ticks_label)

ax_tscale.xaxis.set_minor_formatter(ticker.NullFormatter())

fig_tscale.savefig("Timescale-%s.pdf" % (prefix))



fig_tratio = plt.figure(figsize=(size_fig, size_fig*aratio))
ax_tratio = fig_tratio.add_subplot(111)

color = "k"
ax_tratio.plot(temp_show, texp_med*rcap_med, ls="solid", marker="none", color=color, alpha=0.7, label="Expansion")
ax_tratio.plot(temp_show, avtexp_med*rcap_med, ls="solid", marker="none", color="grey", alpha=0.7, label="Expansion2")

ax_tratio.grid(lw=1.0, ls="dotted", color="grey", alpha=0.5)

ax_tratio.set_xlim(temp_max,temp_min)
ax_tratio.set_xscale("log")
ax_tratio.set_xlabel(temp_label)

ax_tratio.set_ylim(1e-2, 1e2)
ax_tratio.set_yscale("log")
ax_tratio.set_ylabel("$t_\\mathrm{exp}/t_\\mathrm{cap}$")

fig_tratio.subplots_adjust(left=0.15,bottom=0.15,right=0.95, top=0.95)
ax_tratio.legend(fontsize=18,frameon=False,ncols=1,loc='lower right')
ax_tratio.tick_params(axis='both', which='both', top=True, right=True, pad=10)

ax_tratio.set_xticks(temp_ticks_pos)
ax_tratio.set_xticklabels(temp_ticks_label)

ax_tratio.xaxis.set_minor_formatter(ticker.NullFormatter())

fig_tratio.savefig("tratio-%s.pdf" % (prefix))




fig_vr = plt.figure(figsize=(size_fig, size_fig*aratio))
ax_vr = fig_vr.add_subplot(111)

color = "k"
ax_vr.fill_between(temp_show, vr_lw, vr_rw, ls="solid", color=color, alpha=0.2, label="")
ax_vr.plot(temp_show, vr_med, ls="solid", marker="none", color=color, alpha=0.7, label="Expansion")


ax_vr.set_xlim(temp_max,temp_min)
ax_vr.set_xscale("log")
ax_vr.set_xlabel(temp_label)

ax_vr.set_ylim(1e6, 1e10)
ax_vr.set_yscale("log")
ax_vr.set_ylabel("$v^r$ (cm\\,s$^{-1}$)")

fig_vr.subplots_adjust(left=0.15,bottom=0.15,right=0.95, top=0.95)
ax_vr.legend(fontsize=18,frameon=False,ncols=1,loc='lower right')
ax_vr.tick_params(axis='both', which='both', top=True, right=True, pad=10)

ax_vr.set_xticks(temp_ticks_pos)
ax_vr.set_xticklabels(temp_ticks_label)

ax_vr.xaxis.set_minor_formatter(ticker.NullFormatter())

fig_vr.savefig("vr-%s.pdf" % (prefix))



fig_r = plt.figure(figsize=(size_fig, size_fig*aratio))
ax_r = fig_r.add_subplot(111)

color = "k"
ax_r.fill_between(temp_show, r_lw, r_rw, ls="solid", color=color, alpha=0.2, label="")
ax_r.plot(temp_show, r_med, ls="solid", marker="none", color=color, alpha=0.7, label="Expansion")


ax_r.set_xlim(temp_max, temp_min)
ax_r.set_xscale("log")
ax_r.set_xlabel(temp_label)

ax_r.set_xticks(temp_ticks_pos)
ax_r.set_xticklabels(temp_ticks_label)

ax_r.set_ylim(1e6, 3e9)
ax_r.set_yscale("log")
ax_r.set_ylabel("$r$ (cm)")

fig_r.subplots_adjust(left=0.15,bottom=0.15,right=0.95, top=0.95)
ax_r.legend(fontsize=18,frameon=False,ncols=1,loc='lower right')
ax_r.tick_params(axis='both', which='both', top=True, right=True, pad=10)


ax_r.xaxis.set_minor_formatter(ticker.NullFormatter())

fig_r.savefig("r-%s.pdf" % (prefix))


fig_eta = plt.figure(figsize=(size_fig, size_fig*aratio))
ax_eta = fig_eta.add_subplot(111)

color = "k"
ax_eta.fill_between(temp_show, eta_lw, eta_rw, ls="solid", color=color, alpha=0.2, label="")
ax_eta.plot(temp_show, eta_med, ls="solid", marker="none", color=color, alpha=0.7, label="")

ax_eta.set_xlim(temp_max, temp_min)
ax_eta.set_xscale("log")
ax_eta.set_xlabel(temp_label)

ax_eta.set_ylim(0.0,3.0)
ax_eta.set_yscale("linear")
ax_eta.set_ylabel("$\\eta_\\mathrm{e}$")

fig_eta.subplots_adjust(left=0.15,bottom=0.15,right=0.95, top=0.95)
ax_eta.legend(fontsize=18,frameon=False,ncols=1,loc='lower right')
ax_eta.tick_params(axis='both', which='both', top=True, right=True, pad=10)

ax_eta.set_xticks(temp_ticks_pos)
ax_eta.set_xticklabels(temp_ticks_label)

ax_eta.xaxis.set_minor_formatter(ticker.NullFormatter())

#ax_eta.yaxis.set_major_locator(ticker.MultipleLocator(0.1))
#ax_eta.yaxis.set_minor_locator(ticker.MultipleLocator(0.02))


fig_eta.savefig("eta-%s.pdf" % (prefix))


fig_time = plt.figure(figsize=(size_fig, size_fig*aratio))
ax_time = fig_time.add_subplot(111)

color = "k"
ax_time.fill_between(temp_show, time_lw, time_rw, ls="solid", color=color, alpha=0.2, label="")
ax_time.plot(temp_show, time_med, ls="solid", marker="none", color=color, alpha=0.7, label="")

ax_time.set_xlim(temp_max,temp_min)
ax_time.set_xscale("log")
ax_time.set_xlabel(temp_label)

ax_time.set_ylim(0.0,0.3)
ax_time.set_yscale("linear")
ax_time.set_ylabel("$t$ (s)")

fig_time.subplots_adjust(left=0.15,bottom=0.15,right=0.95, top=0.95)
#ax_time.legend(fontsize=18,frameon=False,ncols=1,loc='lower right')
ax_time.tick_params(axis='both', which='both', top=True, right=True, pad=10)

ax_time.set_xticks(temp_ticks_pos)
ax_time.set_xticklabels(temp_ticks_label)

ax_time.xaxis.set_minor_formatter(ticker.NullFormatter())

fig_vr2 = plt.figure(figsize=(size_fig, size_fig*aratio))
ax_vr2 = fig_vr2.add_subplot(111)
fig_temp2 = plt.figure(figsize=(size_fig, size_fig*aratio))
ax_temp2 = fig_temp2.add_subplot(111)

fig_ye_texp = plt.figure(figsize=(size_fig, size_fig*aratio))
ax_ye_texp = fig_ye_texp.add_subplot(111)

str1 = ""
for pid in list_ids:
    if pid in ids_selected:
        str1 += ",%08d" % (pid)
        idx = id_to_index[pid]
        
        dlogT = np.diff(np.log(temp_show))
        dt    = np.diff(time_temp[idx])
        dlogT_dt = dlogT/dt
        dlogT_dt = np.append(dlogT_dt, dlogT_dt[-1])
        
        #dlogT_dt = np.gradient(np.log(temp_show), time_temp[ip])
        #print(1.0/dlogT_dt)
        #ax_tscale.plot(temp_show, 1.0/np.abs(dlogT_dt), ls="solid", lw=0.5, marker="none", color="k", alpha=0.3, label="")

        ax_time.plot(temp_show, time_temp[idx], ls="solid", lw=1.0, marker=".", color="r", alpha=0.4, label="")

        #ax_ye.plot(temp_show, ye_temp[idx], ls="solid", lw=1.0, marker=".", color="k", alpha=0.4, label="")
        #ax_ye.plot(temp_show, yecap_temp[idx], ls="solid", lw=1.0, marker=".", color="r", alpha=0.4, label="")

        #ax_vr.plot(temp_show, vr_temp[idx], ls="solid", lw=1.0, marker=".", color="k", alpha=0.4, label="")
        
        
        fn = "../traj/traj_%08d.dat" % (pid)
        traj_file_exists = os.path.exists(fn)
        print(fn, traj_file_exists)
        if traj_file_exists:
            data = np.loadtxt(fn, comments="#", usecols=[0,1,2,3,4,5,6,8,9])
            t = data[:,0]
            x = data[:,1]
            y = data[:,2]
            z = data[:,3]
            vx= data[:,4]
            vy= data[:,5]
            vz= data[:,6]
            T = data[:,7]
            ye = data[:,8]
            
            r = np.sqrt(x**2 + y**2 + z**2)
            vr= (vx*x + vy*y + vz*z) / r

            time_temp_single = time_temp[idx]
            time_min = np.amin(time_temp_single[ (np.isfinite(time_temp_single)) & (temp_show<temp_max)])

            if show_in_GK:
                T_show = T/1e9
            else:
                T_show = T/1.16e10
                pass
            mask = t > time_min
            ax_time.plot(T_show[mask], t[mask], ls="solid", lw=0.5, marker="none", color="r", alpha=0.4, label="")
            #ax_vr.plot(T[mask]/1e9, vr[mask], ls="solid", lw=0.5, marker="none", color="k", alpha=0.4, label="")
            ax_vr2.plot(t[mask], vr[mask], ls="solid", lw=0.5, marker="none", color="k", alpha=0.4, label="")
            ax_temp2.plot(t[mask], T_show[mask], ls="solid", lw=0.5, marker="none", color="k", alpha=0.4, label="")
            #ax_r.plot(T[mask]/1e9, r[mask], ls="solid", lw=0.5, marker="none", color="k", alpha=0.4, label="")

            ax_ye.plot(T_show[mask], ye[mask], ls="solid", lw=0.5, marker="none", color="k", alpha=0.4, label="")

            ax_ye_texp.plot((r/vr)[mask], ye[mask], ls="solid", lw=0.5, marker="none", color="k", alpha=0.4, label="")
            
            pass

        fn = "../weak/weak_%08d.dat" % (pid)
        print(fn)
        weak_file_exists = os.path.exists(fn)
        if weak_file_exists:
            data = np.loadtxt(fn, comments="#", usecols=[4,5,6,7,8])
            # data = np.loadtxt(
            #     fn,
            #     comments="#",
            #     usecols=[4, 5, 6, 7, 8],
            #     converters=read_fortran_float,
            # )
            eta = data[:,0]
            Xn = data[:,1]
            Xp = data[:,2]
            Rec = data[:,3]
            Rpc = data[:,4]

            # t_weak = 1.0/np.maximum(Rec/Xp,Rpc/Xn)
            t_weak = ye/np.maximum(Rec,Rpc)
            
            if traj_file_exists:
                t_exp  = r/np.abs(vr)
                ax_tratio.plot(T_show[mask], (t_exp/t_weak)[mask], ls="solid", lw=0.5, marker="none", color="k", alpha=0.4, label="")

                ax_tscale.plot(T_show[mask], t_weak[mask], ls="solid", lw=0.5,  marker="none", color="r", alpha=0.4, label="")
                ax_tscale.plot(T_show[mask], t_exp[mask], ls="solid", lw=0.5, marker="none", color="k", alpha=0.4, label="")

                
                pass
            pass
        pass

    pass
print(str1)

for pid in list_ids:
    # if ip in ids_selected:
    idx = id_to_index[pid]
    # for it in range(len(temp_show)):
    #     print(time_temp[idx,it],temp_show[it])
    
    pass
#ax_time.yaxis.set_major_locator(ticker.MultipleLocator(0.1))
#ax_time.yaxis.set_minor_locator(ticker.MultipleLocator(0.02))


fig_time.savefig("time-%s.pdf" % (prefix))

t_min = 0.010; t_max=0.5

ax_vr2.set_xlim(t_min, t_max)
#ax_vr2.set_xscale("log")
ax_vr2.set_xlabel("$t$ (s)")

ax_vr2.set_ylim(-3e9, 3e9)
#ax_vr2.set_yscale("log")
ax_vr2.set_ylabel("$v^r$ (cm\\,s$^{-1}$)")

fig_vr2.subplots_adjust(left=0.15,bottom=0.15,right=0.95, top=0.95)
#ax_vr2.legend(fontsize=18,frameon=False,ncols=1,loc='lower right')
ax_vr2.tick_params(axis='both', which='both', top=True, right=True, pad=10)

#ax_vr2.set_xticks(temp_ticks_pos)
#ax_vr2.set_xticklabels(temp_ticks_label)

ax_vr2.xaxis.set_minor_formatter(ticker.NullFormatter())

fig_vr2.savefig("vr2-%s.pdf" % (prefix))


ax_temp2.set_xlim(t_min, t_max)
#ax_temp2.set_xscale("log")
ax_temp2.set_xlabel("$t$ (s)")

ax_temp2.set_ylim(temp_min, temp_max)
ax_temp2.set_yscale("log")
ax_temp2.set_ylabel(temp_label)
ax_temp2.set_yticks(temp_ticks_pos)
ax_temp2.set_yticklabels(temp_ticks_label)
ax_temp2.yaxis.set_minor_formatter(ticker.NullFormatter())

fig_temp2.subplots_adjust(left=0.15,bottom=0.15,right=0.95, top=0.95)
#ax_temp2.legend(fontsize=18,frameon=False,ncols=1,loc='lower right')
ax_temp2.tick_params(axis='both', which='both', top=True, right=True, pad=10)

fig_temp2.savefig("temp2-%s.pdf" % (prefix))


ax_ye_texp.set_xlim(1e-3,1e3)
ax_ye_texp.set_xscale("log")
ax_ye_texp.set_xlabel("$t_\\mathrm{exp}$ (s)")

ax_ye_texp.set_ylim(0.0,0.6)
#ax_ye_texp.set_yscale("log")
ax_ye_texp.set_ylabel("$Y_\\mathrm{e}$")
#ax_ye_texp.set_yticks(temp_ticks_pos)
#ax_ye_texp.set_yticklabels(temp_ticks_label)
#ax_ye_texp.yaxis.set_minor_formatter(ticker.NullFormatter())

fig_ye_texp.subplots_adjust(left=0.15,bottom=0.15,right=0.95, top=0.95)
#ax_ye_texp.legend(fontsize=18,frameon=False,ncols=1,loc='lower right')
ax_ye_texp.tick_params(axis='both', which='both', top=True, right=True, pad=10)

fig_ye_texp.savefig("ye-t_exp-%s.pdf" % (prefix))


plt.show()

