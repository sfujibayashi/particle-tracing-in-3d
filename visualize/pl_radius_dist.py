import numpy as np
import h5py
import sys
#import subprocess
import os
import glob
import argparse

from scipy.interpolate import interp2d

import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as colors
import matplotlib.ticker as ticker
from matplotlib.colors import LinearSegmentedColormap

from mpl_toolkits import mplot3d

# plt.switch_backend('agg')
# os.environ["PATH"] += os.pathsep + '/data/home/sfujibayashi/libs/texlive/2020/bin/x86_64-linux'

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

msun = 1.988e33
clight=2.99792458e10
MeV_to_K = 1.16e10

parser = argparse.ArgumentParser()
parser.add_argument('--dir_read', required=True, type=str)

args = parser.parse_args()
dir_read = args.dir_read

list_file = sorted(glob.glob("%s/data_???_??????.h5" % (dir_read)))
# print(list_file)

# list_time = []
time_1 = 1e99
time_2 = 0.0
for fn_data in list_file:
    f1 = h5py.File(fn_data, 'r')
    time = np.asarray(f1['/time'])[0]
    
    time_1 = np.min([time_1, time])
    time_2 = np.max([time_2, time])
    pass

cmap_time = plt.get_cmap('jet', 256)
norm_time = colors.Normalize(vmin=time_1, vmax=time_2)


fig_radius = plt.figure(figsize=(10.0, 10.0*0.624))
ax_radius = fig_radius.add_subplot(111)

fig_temp = plt.figure(figsize=(10.0, 10.0*0.624))
ax_temp = fig_temp.add_subplot(111)

bins_r = np.logspace(7.0, 9.0, 101)
bins_T = np.logspace(9.0, 11.0, 101)

alpha=0.7

T_thr = 7.0e9; str1="7GK"
R_thr = 3.0e8


for fn_data in list_file:
    f1 = h5py.File(fn_data, 'r')

    ntraj = np.asarray(f1['/np'])[0]
    time = np.asarray(f1['/time'])[0]

    print("# of tracers = %10d, time = %e" % (ntraj,time))
    
    dm_p = np.asarray(f1['/dm_p'])/msun
    x_p = np.asarray(f1['/x_p'])
    y_p = np.asarray(f1['/y_p'])
    z_p = np.asarray(f1['/z_p'])
    vlx_p = np.asarray(f1['/vlx_p'])
    vly_p = np.asarray(f1['/vly_p'])
    vlz_p = np.asarray(f1['/vlz_p'])
    
    tem_p = np.asarray(f1['/tem_p'])*MeV_to_K
    ut1_p = np.asarray(f1['/ut1_p'])
    ye_p = np.asarray(f1['/ye_p'])
    entr_p = np.asarray(f1['/sen_p'])
    
    r_p   = np.sqrt(x_p**2 + y_p**2 + z_p**2)
    # vlr_p = np.sqrt(vlx_p**2 + vly_p**2 + vlz_p**2)*clight
    vlr_p = (vlx_p*x_p + vly_p*y_p + vlz_p*z_p)/r_p * clight
    
    mass_lt3e8cm = np.sum(dm_p[r_p < 3.0e8])
    mass_gt3GK   = np.sum(dm_p[tem_p > 3.0e9])
    label_r = "$M(r<3\\times10^8\\mathrm{cm})=%6.3f M_\\odot$" % (mass_lt3e8cm)
    label_T = "$M(T>3\\mathrm{GK})=%6.3f M_\\odot$" % (mass_gt3GK)
    ####

    color = cmap_time(norm_time(time))

    counts, edges = np.histogram(r_p, weights=dm_p, bins=bins_r)
    cdf = np.cumsum(counts)
    ax_radius.stairs(cdf, edges, baseline=None, linewidth=3.0, color=color, label=label_r)


    counts, edges = np.histogram(tem_p, weights=dm_p, bins=bins_T)
    cdf = np.cumsum(counts[::-1])[::-1]
    ax_temp.stairs(cdf, edges, baseline=None, linewidth=3.0, color=color, label=label_T)

                     
    #h = ax_radius.hist(r_p, weights=dm_p,
    #                   bins=bins, log=False, alpha=alpha, align='mid', facecolor="w", edgecolor="grey", histtype="step", cumulative=True, lw=3.0, label="$hu_t+h_\\mathrm{min}<0$ ($%10.5fM_\\odot$)" % (mass_tot))
    pass

ax_radius.grid(ls="dotted", color = "grey", alpha=0.5)

s_map = cm.ScalarMappable(norm=norm_time, cmap=cmap_time)
s_map.set_array([])
cbar=fig_radius.colorbar(s_map,pad=0.01, ax=ax_radius)
#cbar.set_ticks([6.0/6.0, 7.0/6.0, 8.0/6.0, 9.0/6.0, 10.0/6.0])
#cbar.ax.set_yticklabels(["1", "7/6", "4/3", "3/2", "5/3"])
cbar.set_label("$t$ (s)",rotation = -90, labelpad= 30)


fig_radius.subplots_adjust(left=0.15,bottom=0.15,right=0.95, top=0.95)
ax_radius.legend(fontsize=12, frameon=True, ncols=1, loc='upper left')
ax_radius.tick_params(axis='both', which='both', top=True, right=True)
ax_radius.set_xlabel("$r$ (cm)")
ax_radius.set_xscale('log')
#ax_radius.set_xlim(3.0,3e2)
    
#ax_radius.xaxis.set_major_locator(ticker.MultipleLocator(0.1))
#ax_radius.xaxis.set_minor_locator(ticker.MultipleLocator(0.01))

ax_radius.set_ylabel("$M(<r)/M_\\odot$")
#ax_radius.set_ylim(1e-5,1.0)
#ax_radius.set_yscale('log')
#ax_radius.yaxis.set_major_locator(ticker.LogLocator(numticks=20))
#ax_radius.yaxis.set_minor_locator(ticker.LogLocator(numticks=20, subs=(.1,.2,.3,.4,.5,.6,.7,.8,.9)))

fig_radius.savefig("radius_hist.pdf")



ax_temp.grid(ls="dotted", color = "grey", alpha=0.5)

s_map = cm.ScalarMappable(norm=norm_time, cmap=cmap_time)
s_map.set_array([])
cbar=fig_temp.colorbar(s_map,pad=0.01, ax=ax_temp)
#cbar.set_ticks([6.0/6.0, 7.0/6.0, 8.0/6.0, 9.0/6.0, 10.0/6.0])
#cbar.ax.set_yticklabels(["1", "7/6", "4/3", "3/2", "5/3"])
cbar.set_label("$t$ (s)",rotation = -90, labelpad= 30)


fig_temp.subplots_adjust(left=0.15,bottom=0.15,right=0.95, top=0.95)
ax_temp.legend(fontsize=12, frameon=True, ncols=1, loc='lower right')
ax_temp.tick_params(axis='both', which='both', top=True, right=True)
ax_temp.set_xlabel("$T$ (K)")
ax_temp.set_xscale('log')
#ax_temp.set_xlim(3.0,3e2)
    
#ax_temp.xaxis.set_major_locator(ticker.MultipleLocator(0.1))
#ax_temp.xaxis.set_minor_locator(ticker.MultipleLocator(0.01))

ax_temp.set_ylabel("$M(>T)/M_\\odot$")
#ax_temp.set_ylim(1e-5,1.0)
#ax_temp.set_yscale('log')
#ax_temp.yaxis.set_major_locator(ticker.LogLocator(numticks=20))
#ax_temp.yaxis.set_minor_locator(ticker.LogLocator(numticks=20, subs=(.1,.2,.3,.4,.5,.6,.7,.8,.9)))

fig_temp.savefig("temp_hist.pdf")


plt.show()
