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

grav = 6.674e-8
msun = 1.988e33
clight=2.99792458e10
MeV_to_K = 1.16e10

parser = argparse.ArgumentParser()
parser.add_argument('--file', required=True, type=str)
parser.add_argument('--file_ref', required=False, type=str, default = "")


args = parser.parse_args()
fn_data = args.file
fn_ref = args.file_ref

i = fn_data.rfind("/")
file_name = fn_data[i+1:]
if file_name[0:5]=="data":
    job= int(file_name[5:8])
    it = int(file_name[9:15])
    label_file = "%03d_%06d" % (job,it)
else:
    job= int(file_name[5:7])
    label_file = "%03d_last" % (job)
    pass

T_thr = 7.0e9; str1="7GK"
R_thr = 3.0e8

    
f1 = h5py.File(fn_data, 'r')

ntraj = np.asarray(f1['/np'])[0]
time = np.asarray(f1['/time'])[0]

print(ntraj,time)

dm_p = np.asarray(f1['/dm_p'])/msun
x_p = np.asarray(f1['/x_p'])
y_p = np.asarray(f1['/y_p'])
z_p = np.asarray(f1['/z_p'])
vlx_p = np.asarray(f1['/vlx_p'])
vly_p = np.asarray(f1['/vly_p'])
vlz_p = np.asarray(f1['/vlz_p'])

tem_p = np.asarray(f1['/tem_p'])*MeV_to_K
ut1_p = np.asarray(f1['/ut1_p'])
hut_p = np.asarray(f1['/hut_p'])
ut_p = np.asarray(f1['/ut_p'])
hhh_p = np.asarray(f1['/hhh_p'])
ye_p = np.asarray(f1['/ye_p'])
entr_p = np.asarray(f1['/sen_p'])

rcyl_p = np.sqrt(x_p**2 + y_p**2)
r_p   = np.sqrt(x_p**2 + y_p**2 + z_p**2)
polar_angle_p = np.degrees(np.arctan2(rcyl_p, z_p))

vlr_p = (vlx_p*x_p + vly_p*y_p + vlz_p*z_p)/r_p * clight
# vlr_p = np.maximum(vlr_p, 1.1e6)
vphi_p = (-vlx_p*y_p + vly_p*x_p)/np.sqrt(x_p**2+y_p**2) * clight

hhh_min = 0.998586047091485
hut_now = hhh_p*ut_p + hhh_min
ut1_now = ut_p + 1.0


mass_tot = np.sum(dm_p)
mass_over_thr = np.sum(dm_p[tem_p >= T_thr])
mass_under_thr = np.sum(dm_p[tem_p < T_thr])
mass_under_thr_ut1 = np.sum(dm_p[(tem_p < T_thr) & (ut1_now<0.0)])
mass_ut1 = np.sum(dm_p[ut1_now<0.0])
mass_hut = np.sum(dm_p[hut_now<0.0])

mass_under_R_thr = np.sum(dm_p[r_p<R_thr])
mass_bet_R = np.sum(dm_p[(3.0e8 < r_p) & (r_p<4.0e8)])
mass_bet_R2 = np.sum(dm_p[(2.0e8 < r_p) & (r_p<3.0e8)])

print("T_thr             =",T_thr)
print("M(tot)            =",mass_tot)
print("M(T>=thr)         =",mass_over_thr)
print("M(T< thr)         =",mass_under_thr)
print("M(T< thr & ut+1<0)=",mass_under_thr_ut1)
print("M(R< thr)=",mass_under_R_thr)

mass_become_ejecta_bernoulli = np.sum(dm_p[(hut_p > 0.0) & (hut_now<0.0)])
print("M(become ejecta)  =",mass_become_ejecta_bernoulli)


mass_over_R_thr = np.sum(dm_p[r_p > R_thr])
print("M(become over R_thr)  =", mass_over_R_thr)



list_ecolors = ["r" if ut1_now[i]<0.0 else "k" for i in range(ntraj) ]

fig_r = plt.figure(figsize=(10.0, 10.0*0.70))
ax_r = fig_r.add_subplot(111)

skip_plot = 1
ax_r.scatter(r_p[::skip_plot*10], tem_p[::skip_plot*10], marker=".", color="grey", edgecolor="k", label="All ($%10.5fM_\\odot$)" % (mass_tot), zorder=0.0)
ax_r.scatter(r_p[hut_now<0.0][::skip_plot], tem_p[hut_now<0.0][::skip_plot], marker="o", color="orange", edgecolor="k", label="$hu_t+h_\\mathrm{min}<0$ ($%10.5fM_\\odot$)" % (mass_hut), zorder=1.0)
#ax_r.scatter(r_p[hut_now<0.0][::skip_plot], tem_p[hut_now<0.0][::skip_plot], marker="o", edgecolor=list_ecolors[::skip_plot], label="$u_t+1<0$ ($%10.5fM_\\odot$)" % (mass_ut1), zorder=1.0)
ax_r.scatter(r_p[ut1_now<0.0][::skip_plot], tem_p[ut1_now<0.0][::skip_plot], marker="o", color="orange", edgecolor="r", label="$u_t+1<0$ ($%10.5fM_\\odot$)" % (mass_ut1), zorder=2.0)

ax_r.text(0.02, 0.10, "$M(T\\geq\\mathrm{%s})=%10.5fM_\\odot$" % (str1,mass_over_thr), ha="left", va="bottom", transform=ax_r.transAxes, fontsize=14)
ax_r.text(0.02, 0.06, "$M(T<\\mathrm{%s})=%10.5fM_\\odot$" % (str1,mass_under_thr), ha="left", va="bottom", transform=ax_r.transAxes, fontsize=14)
ax_r.text(0.02, 0.02, "$M(T<\\mathrm{%s}\&u_t+1>0)=%10.5fM_\\odot$" % (str1,mass_under_thr-mass_under_thr_ut1), ha="left", va="bottom", transform=ax_r.transAxes, fontsize=14)

# ax_r.plot([1e6,5e8], [T_thr]*2, ls="dotted", color="grey", lw=1.0)

ax_r.grid(ls="dotted", color="grey", lw=1.0)

ax_r.text(0.98, 0.02, s="$t=%6.4f$ s" % (time), ha="right", va="bottom", transform=ax_r.transAxes)

fig_r.subplots_adjust(left=0.15,bottom=0.15,right=0.95, top=0.95)
ax_r.legend(fontsize=14, frameon=True, ncols=1, loc='upper right')
ax_r.tick_params(axis='both', which='both', top=True, right=True)
ax_r.set_xlabel("$r$ (cm)")
ax_r.set_xscale('log')
ax_r.set_xlim(0.7e6,1e10)
#ax_r.xaxis.set_major_locator(ticker.MultipleLocator(1000))
#ax_r.xaxis.set_minor_locator(ticker.MultipleLocator(100))

ax_r.set_ylabel("$T$ (K)")
ax_r.set_ylim(3e8,1e11)
ax_r.set_yscale('log')
ax_r.yaxis.set_major_locator(ticker.LogLocator(numticks=20))
ax_r.yaxis.set_minor_locator(ticker.LogLocator(numticks=20, subs=(.1,.2,.3,.4,.5,.6,.7,.8,.9)))

fig_r.savefig("r_temp_%s.png" % (label_file))



###

fig_vel = plt.figure(figsize=(10.0, 10.0*0.70))
ax_vel = fig_vel.add_subplot(111)

skip_plot = 1
ax_vel.scatter(r_p[::skip_plot*10], vlr_p[::skip_plot*10], marker=".", color="grey", edgecolor="k", label="All ($%10.5fM_\\odot$)" % (mass_tot), zorder=0.0)
ax_vel.scatter(r_p[hut_now<0.0][::skip_plot], vlr_p[hut_now<0.0][::skip_plot], marker="o", color="orange", edgecolor="k", label="$hu_t+h_\\mathrm{min}<0$ ($%10.5fM_\\odot$)" % (mass_hut), zorder=1.0)
# ax_vel.scatter(r_p[ (hut_now<0.0) & (vlr_p<0.0)][::skip_plot], -vlr_p[(hut_now<0.0)&(vlr_p<0.0)][::skip_plot], marker="o", s=10, color="orange", edgecolor="b", label="$hu_t+h_\\mathrm{min}<0$ \\& $v^r<0$", zorder=1.0)
ax_vel.scatter(r_p[ut1_now<0.0][::skip_plot], vlr_p[ut1_now<0.0][::skip_plot], marker="o", color="orange", edgecolor="r", label="$u_t+1<0$ ($%10.5fM_\\odot$)" % (mass_ut1), zorder=2.0)
# ax_vel.scatter(r_p[(ut1_now<0.0) & (polar_angle_p < 15)][::skip_plot], vlr_p[(ut1_now<0.0) & (polar_angle_p < 15)][::skip_plot], marker="o", color="r", zorder=5.0)

# ax_vel.text(0.02, 0.10, "$M(T\\geq\\mathrm{%s})=%10.5fM_\\odot$" % (mass_over_thr), ha="left", va="bottom", transform=ax_vel.transAxes, fontsize=14)
# ax_vel.text(0.02, 0.06, "$M(T<\\mathrm{%s})=%10.5fM_\\odot$" % (mass_under_thr), ha="left", va="bottom", transform=ax_vel.transAxes, fontsize=14)
# ax_vel.text(0.02, 0.02, "$M(T<\\mathrm{%s}\&u_t+1>0)=%10.5fM_\\odot$" % (mass_under_thr-mass_under_thr_ut1), ha="left", va="bottom", transform=ax_vel.transAxes, fontsize=14)

# ax_vel.plot([1e6,5e8], [T_thr]*2, ls="dotted", color="grey", lw=1.0)

r  = np.array([1e7,1e9])
M0 = 2.5*msun
vkep=np.sqrt(grav*M0/r)
ax_vel.plot(r, vkep, ls="dashed", color="grey", lw=2.0, label="$v^r = v_\\mathrm{K}$ (for $M=2.5M_\\odot$)")

ax_vel.grid(ls="dotted", color="grey", lw=1.0)

ax_vel.text(0.98, 0.02, s="$t=%6.4f$ s" % (time), ha="right", va="bottom", transform=ax_vel.transAxes)


fig_vel.subplots_adjust(left=0.15,bottom=0.15,right=0.95, top=0.95)
ax_vel.legend(fontsize=14, frameon=True, ncols=1, loc='lower left')
ax_vel.tick_params(axis='both', which='both', top=True, right=True)
ax_vel.set_xlabel("$r$ (cm)")
ax_vel.set_xscale('log')
ax_vel.set_xlim(0.7e6,1e10)
#ax_vel.xaxis.set_major_locator(ticker.MultipleLocator(1000))
#ax_vel.xaxis.set_minor_locator(ticker.MultipleLocator(100))

ax_vel.set_ylabel("$v^r$ (cm/s)")
ax_vel.set_ylim(1e6,3e10)
ax_vel.set_yscale('log')
ax_vel.yaxis.set_major_locator(ticker.LogLocator(numticks=20))
ax_vel.yaxis.set_minor_locator(ticker.LogLocator(numticks=20, subs=(.1,.2,.3,.4,.5,.6,.7,.8,.9)))

fig_vel.savefig("r_vel_%s.png" % (label_file))


###

fig_vphi = plt.figure(figsize=(10.0, 10.0*0.70))
ax_vphi = fig_vphi.add_subplot(111)

skip_plot = 1
ax_vphi.scatter(r_p[::skip_plot*10], vphi_p[::skip_plot*10], marker=".", color="grey", edgecolor="k", label="All ($%10.5fM_\\odot$)" % (mass_tot), zorder=0.0)
ax_vphi.scatter(r_p[hut_now<0.0][::skip_plot], vphi_p[hut_now<0.0][::skip_plot], marker="o", color="orange", edgecolor="k", label="$hu_t+h_\\mathrm{min}<0$ ($%10.5fM_\\odot$)" % (mass_hut), zorder=1.0) 
ax_vphi.scatter(r_p[ut1_now<0.0][::skip_plot], vphi_p[ut1_now<0.0][::skip_plot], marker="o", color="orange", edgecolor="r", label="$u_t+1<0$ ($%10.5fM_\\odot$)" % (mass_ut1), zorder=2.0)

# ax_vphi.text(0.02, 0.10, "$M(T\\geq\\mathrm{%s})=%10.5fM_\\odot$" % (mass_over_thr), ha="left", va="bottom", transform=ax_vphi.transAxes, fontsize=14)
# ax_vphi.text(0.02, 0.06, "$M(T<\\mathrm{%s})=%10.5fM_\\odot$" % (mass_under_thr), ha="left", va="bottom", transform=ax_vphi.transAxes, fontsize=14)
# ax_vphi.text(0.02, 0.02, "$M(T<\\mathrm{%s}\&u_t+1>0)=%10.5fM_\\odot$" % (mass_under_thr-mass_under_thr_ut1), ha="left", va="bottom", transform=ax_vphi.transAxes, fontsize=14)

# ax_vphi.plot([1e6,5e8], [T_thr]*2, ls="dotted", color="grey", lw=1.0)

r  = np.array([1e7,1e9])
M0 = 2.5*msun
vkep=np.sqrt(grav*M0/r)
ax_vphi.plot(r, vkep, ls="dotted", color="grey", lw=1.0, label="$v^r = v_\\mathrm{K}$ (for $M=2.5M_\\odot$)")

fig_vphi.subplots_adjust(left=0.15,bottom=0.15,right=0.95, top=0.95)
ax_vphi.legend(fontsize=14, frameon=True, ncols=1, loc='lower left')
ax_vphi.tick_params(axis='both', which='both', top=True, right=True)
ax_vphi.set_xlabel("$r$ (cm)")
ax_vphi.set_xscale('log')
ax_vphi.set_xlim(0.7e6,1e10)
#ax_vphi.xaxis.set_major_locator(ticker.MultipleLocator(1000))
#ax_vphi.xaxis.set_minor_locator(ticker.MultipleLocator(100))

ax_vphi.set_ylabel("$v^\\phi$ (cm/s)")
ax_vphi.set_ylim(1e6,3e10)
ax_vphi.set_yscale('log')
ax_vphi.yaxis.set_major_locator(ticker.LogLocator(numticks=20))
ax_vphi.yaxis.set_minor_locator(ticker.LogLocator(numticks=20, subs=(.1,.2,.3,.4,.5,.6,.7,.8,.9)))

fig_vphi.savefig("r_vphi_%s.png" % (label_file))


###

fig_vratio = plt.figure(figsize=(10.0, 10.0*0.70))
ax_vratio = fig_vratio.add_subplot(111)

skip_plot = 1
ax_vratio.scatter(r_p[::skip_plot*10], (vlr_p/vphi_p)[::skip_plot*10], marker=".", color="grey", edgecolor="k", label="All ($%10.5fM_\\odot$)" % (mass_tot), zorder=0.0)
ax_vratio.scatter(r_p[hut_now<0.0][::skip_plot], (vlr_p/vphi_p)[hut_now<0.0][::skip_plot], marker="o", color="orange", edgecolor="k", label="$hu_t+h_\\mathrm{min}<0$ ($%10.5fM_\\odot$)" % (mass_hut), zorder=1.0)
ax_vratio.scatter(r_p[ut1_now<0.0][::skip_plot], (vlr_p/vphi_p)[ut1_now<0.0][::skip_plot], marker="o", color="orange", edgecolor="r", label="$u_t+1<0$ ($%10.5fM_\\odot$)" % (mass_ut1), zorder=2.0) 

# ax_vratio.text(0.02, 0.10, "$M(T\\geq\\mathrm{%s})=%10.5fM_\\odot$" % (mass_over_thr), ha="left", va="bottom", transform=ax_vratio.transAxes, fontsize=14)
# ax_vratio.text(0.02, 0.06, "$M(T<\\mathrm{%s})=%10.5fM_\\odot$" % (mass_under_thr), ha="left", va="bottom", transform=ax_vratio.transAxes, fontsize=14)
# ax_vratio.text(0.02, 0.02, "$M(T<\\mathrm{%s}\&u_t+1>0)=%10.5fM_\\odot$" % (mass_under_thr-mass_under_thr_ut1), ha="left", va="bottom", transform=ax_vratio.transAxes, fontsize=14)

# ax_vratio.plot([1e6,5e8], [T_thr]*2, ls="dotted", color="grey", lw=1.0)
ax_vratio.grid(ls="dotted", color = "grey", alpha=0.5)

fig_vratio.subplots_adjust(left=0.15,bottom=0.15,right=0.95, top=0.95)
ax_vratio.legend(fontsize=14, frameon=True, ncols=1, loc='upper left')
ax_vratio.tick_params(axis='both', which='both', top=True, right=True)
ax_vratio.set_xlabel("$r$ (cm)")
ax_vratio.set_xscale('log')
ax_vratio.set_xlim(0.7e6,1e10)
#ax_vratio.xaxis.set_major_locator(ticker.MultipleLocator(1000))
#ax_vratio.xaxis.set_minor_locator(ticker.MultipleLocator(100))

ax_vratio.set_ylabel("$v^r/v^\\phi$")
ax_vratio.set_ylim(0.01,100.0)
ax_vratio.set_yscale('log')
ax_vratio.yaxis.set_major_locator(ticker.LogLocator(numticks=20))
ax_vratio.yaxis.set_minor_locator(ticker.LogLocator(numticks=20, subs=(.1,.2,.3,.4,.5,.6,.7,.8,.9)))

fig_vratio.savefig("r_vratio_%s.png" % (label_file))


###

fig_ye = plt.figure(figsize=(10.0, 10.0*0.624))
ax_ye = fig_ye.add_subplot(111)

bins = np.linspace(0.0, 0.6, 61)
#print(bins)
#sys.exit()
alpha=0.7

ax_ye.hist(ye_p, weights=dm_p/mass_tot,
           bins=bins, log=True, alpha=alpha, align='mid', facecolor="w", edgecolor="grey", histtype="step", lw=3.0, label="$hu_t+h_\\mathrm{min}<0$ ($%10.5fM_\\odot$)" % (mass_tot))
# ax_ye.hist(ye_p[tem_p < T_thr], weights=dm_p[tem_p < T_thr]/mass_tot,
#            bins=bins, log=True, alpha=alpha, align='mid', facecolor="w", edgecolor="r", histtype="step", lw=1.5, label="$T<\\mathrm{%s}$ ($%10.5fM_\\odot$)" % (str1,mass_under_thr))
# ax_ye.hist(ye_p[tem_p >= T_thr], weights=dm_p[tem_p >= T_thr]/mass_tot,
#            bins=bins, log=True, alpha=alpha, align='mid', facecolor="w", edgecolor="b", histtype="step", lw=1.5, label="$T\\geq\\mathrm{%s}$ ($%10.5fM_\\odot$)" % (str1,mass_over_thr))
# ax_ye.hist(ye_p[(tem_p < T_thr) & (ut1_now>0.0)], weights=dm_p[(tem_p < T_thr) & (ut1_now>0.0)]/mass_tot,
#            bins=bins, log=True, alpha=alpha, align='mid', facecolor="w", edgecolor="g", histtype="step", lw=1.5, label="$T<\\mathrm{%s}\&u_t+1>0$ ($%10.5fM_\\odot$)" % (str1,mass_under_thr-mass_under_thr_ut1))
# ax_ye.hist(ye_p[ut1_now<0.0], weights=dm_p[ut1_now<0.0]/mass_tot,
#            bins=bins, log=True, alpha=alpha, align='mid', facecolor="w", edgecolor="orange", histtype="step", lw=1.5, ls="dashed", label="$u_t+1<0$ ($%10.5fM_\\odot$)" % (mass_ut1))
ax_ye.hist(ye_p[r_p<3.0e8], weights=dm_p[r_p<3.0e8]/mass_tot,
           bins=bins, log=True, alpha=alpha, align='mid', facecolor="w", edgecolor="cyan", histtype="step", lw=1.5, ls="solid", label="$r<3\\times 10^8$cm ($%10.5fM_\\odot$)" % (mass_under_R_thr))
ax_ye.hist(ye_p[(3.0e8 < r_p) & (r_p<4.0e8)], weights=dm_p[(3.0e8 < r_p) & (r_p<4.0e8)]/mass_tot,
           bins=bins, log=True, alpha=alpha, align='mid', facecolor="w", edgecolor="blue", histtype="step", lw=1.5, ls="solid", label="$3\\times 10^8 \\mathrm{cm} - 4\\times 10^8$cm ($%10.5fM_\\odot$)" % (mass_bet_R))
ax_ye.hist(ye_p[(2.0e8 < r_p) & (r_p<3.0e8)], weights=dm_p[(2.0e8 < r_p) & (r_p<3.0e8)]/mass_tot,
           bins=bins, log=True, alpha=alpha, align='mid', facecolor="w", edgecolor="red", histtype="step", lw=1.5, ls="solid", label="$2\\times 10^8 \\mathrm{cm} - 3\\times 10^8$cm ($%10.5fM_\\odot$)" % (mass_bet_R2))

fig_ye.subplots_adjust(left=0.15,bottom=0.15,right=0.95, top=0.95)
ax_ye.legend(fontsize=14, frameon=True, ncols=1, loc='upper right')
ax_ye.tick_params(axis='both', which='both', top=True, right=True)
ax_ye.set_xlabel("$Y_\\mathrm{e}$")
#ax_ye.set_xscale('log')
ax_ye.set_xlim(0.15,0.6)

ax_ye.xaxis.set_major_locator(ticker.MultipleLocator(0.1))
ax_ye.xaxis.set_minor_locator(ticker.MultipleLocator(0.01))

ax_ye.set_ylabel("$\\Delta M/M_\\mathrm{ej}$")
ax_ye.set_ylim(1e-5,1.0)
ax_ye.set_yscale('log')
ax_ye.yaxis.set_major_locator(ticker.LogLocator(numticks=20))
ax_ye.yaxis.set_minor_locator(ticker.LogLocator(numticks=20, subs=(.1,.2,.3,.4,.5,.6,.7,.8,.9)))

fig_ye.savefig("Ye.pdf")


####

fig_entr = plt.figure(figsize=(10.0, 10.0*0.624))
ax_entr = fig_entr.add_subplot(111)

bins = np.logspace(0.0, 3.0, 61)
#print(bins)
#sys.exit()
alpha=0.7

ax_entr.hist(entr_p, weights=dm_p/mass_tot,
           bins=bins, log=True, alpha=alpha, align='mid', facecolor="w", edgecolor="grey", histtype="step", lw=3.0, label="$hu_t+h_\\mathrm{min}<0$ ($%10.5fM_\\odot$)" % (mass_tot))
ax_entr.hist(entr_p[tem_p < T_thr], weights=dm_p[tem_p < T_thr]/mass_tot,
           bins=bins, log=True, alpha=alpha, align='mid', facecolor="w", edgecolor="r", histtype="step", lw=1.5, label="$T<\\mathrm{%s}$ ($%10.5fM_\\odot$)" % (str1,mass_under_thr))
ax_entr.hist(entr_p[tem_p >= T_thr], weights=dm_p[tem_p >= T_thr]/mass_tot,
           bins=bins, log=True, alpha=alpha, align='mid', facecolor="w", edgecolor="b", histtype="step", lw=1.5, label="$T\\geq\\mathrm{%s}$ ($%10.5fM_\\odot$)" % (str1,mass_over_thr))
ax_entr.hist(entr_p[(tem_p < T_thr) & (ut1_now>0.0)], weights=dm_p[(tem_p < T_thr) & (ut1_now>0.0)]/mass_tot,
           bins=bins, log=True, alpha=alpha, align='mid', facecolor="w", edgecolor="g", histtype="step", lw=1.5, label="$T<\\mathrm{%s}\&u_t+1>0$ ($%10.5fM_\\odot$)" % (str1,mass_under_thr-mass_under_thr_ut1))
ax_entr.hist(entr_p[ut1_now<0.0], weights=dm_p[ut1_now<0.0]/mass_tot,
           bins=bins, log=True, alpha=alpha, align='mid', facecolor="w", edgecolor="orange", histtype="step", lw=1.5, ls="dashed", label="$u_t+1<0$ ($%10.5fM_\\odot$)" % (mass_ut1))

fig_entr.subplots_adjust(left=0.15,bottom=0.15,right=0.95, top=0.95)
ax_entr.legend(fontsize=14, frameon=True, ncols=1, loc='upper right')
ax_entr.tick_params(axis='both', which='both', top=True, right=True)
ax_entr.set_xlabel("$s/k$")
ax_entr.set_xscale('log')
ax_entr.set_xlim(3.0,3e2)

#ax_entr.xaxis.set_major_locator(ticker.MultipleLocator(0.1))
#ax_entr.xaxis.set_minor_locator(ticker.MultipleLocator(0.01))

ax_entr.set_ylabel("$\\Delta M/M_\\mathrm{ej}$")
ax_entr.set_ylim(1e-5,1.0)
ax_entr.set_yscale('log')
ax_entr.yaxis.set_major_locator(ticker.LogLocator(numticks=20))
ax_entr.yaxis.set_minor_locator(ticker.LogLocator(numticks=20, subs=(.1,.2,.3,.4,.5,.6,.7,.8,.9)))

fig_entr.savefig("entr.pdf")


####

fig_temp = plt.figure(figsize=(10.0, 10.0*0.624))
ax_temp = fig_temp.add_subplot(111)

bins = np.logspace(9.0, 11.0, 101)
#print(bins)
#sys.exit()
alpha=0.7

counts, edges = np.histogram(tem_p, weights=dm_p, bins=bins)
cdf = np.cumsum(counts)
ax_temp.stairs(cdf, edges, baseline=None, linewidth=3.0,
               label=r"$hu_t+h_{\mathrm{min}}<0$ ($%10.5fM_\odot$)" % mass_tot)

#ax_temp.hist(tem_p, weights=dm_p,
#             bins=bins, log=False, alpha=alpha, align='mid', facecolor="w", edgecolor="grey", histtype="step", cumulative=True, lw=3.0, label="$hu_t+h_\\mathrm{min}<0$ ($%10.5fM_\\odot$)" % (mass_tot))
ax_temp.grid(ls="dotted", color = "grey", alpha=0.5)

fig_temp.subplots_adjust(left=0.15,bottom=0.15,right=0.95, top=0.95)
ax_temp.legend(fontsize=14, frameon=True, ncols=1, loc='lower right')
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

fig_temp.savefig("temp.pdf")




####

fig_radius = plt.figure(figsize=(10.0, 10.0*0.624))
ax_radius = fig_radius.add_subplot(111)

bins = np.logspace(7.0, 9.0, 101)
#print(bins)
#sys.exit()
alpha=0.7

counts, edges = np.histogram(r_p, weights=dm_p, bins=bins)
cdf = np.cumsum(counts)
ax_radius.stairs(cdf, edges, baseline=None, linewidth=3.0,
                 label=r"$hu_t+h_{\mathrm{min}}<0$ ($%10.5fM_\odot$)" % mass_tot)
#h = ax_radius.hist(r_p, weights=dm_p,
#                   bins=bins, log=False, alpha=alpha, align='mid', facecolor="w", edgecolor="grey", histtype="step", cumulative=True, lw=3.0, label="$hu_t+h_\\mathrm{min}<0$ ($%10.5fM_\\odot$)" % (mass_tot))
ax_radius.grid(ls="dotted", color = "grey", alpha=0.5)

fig_radius.subplots_adjust(left=0.15,bottom=0.15,right=0.95, top=0.95)
ax_radius.legend(fontsize=14, frameon=True, ncols=1, loc='upper left')
ax_radius.tick_params(axis='both', which='both', top=True, right=True)
ax_radius.set_xlabel("$r$ (cm)")
ax_radius.set_xscale('log')
#ax_radius.set_xlim(3.0,3e2)

#ax_radius.xaxis.set_major_locator(ticker.MultipleLocator(0.1))
#ax_radius.xaxis.set_minor_locator(ticker.MultipleLocator(0.01))

ax_radius.set_ylabel("$M(>r)/M_\\odot$")
#ax_radius.set_ylim(1e-5,1.0)
#ax_radius.set_yscale('log')
#ax_radius.yaxis.set_major_locator(ticker.LogLocator(numticks=20))
#ax_radius.yaxis.set_minor_locator(ticker.LogLocator(numticks=20, subs=(.1,.2,.3,.4,.5,.6,.7,.8,.9)))

fig_radius.savefig("radius_hist.pdf")




plt.show()
