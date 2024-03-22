import sys
#sys.path.insert(0, '/sakura/ptmp/lcmenegazzi/athena-public-version/vis/python')
#import athena_read
import numpy as np
import os
import h5py
import argparse
import sys
import glob
import re

import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as colors
import matplotlib.ticker as ticker

def natural_sort(l): 
    convert = lambda text: int(text) if text.isdigit() else text.lower()
    alphanum_key = lambda key: [convert(c) for c in re.split('([0-9]+)', key)]
    return sorted(l, key=alphanum_key)

# solar_radius=6.960012186e10
msun=1.988e33
# t_uni = solar_radius/(1.e8)
# rho_uni = solar_mass/solar_radius**3
# vel_uni = 1e8
# ene_uni = rho_uni*vel_uni**2
size_fig = 10.0

r_max = 3e7

fn = "/sakura/ptmp/shofu/DD2Tim626_135_135_44km_150mstg_B0_HLLC_FUKA/Analysis_ptr/data/amr_steps_6.h5"

dir_out = "/sakura/ptmp/shofu/DD2Tim626_135_135_44km_150mstg_B0_HLLC_FUKA/Analysis_ptr/fig"


f1 = h5py.File(fn, 'r')

list_levels = list(f1["/"].keys())
list_levels = natural_sort(list_levels)
# sys.exit()
list_steps = []
for level in list_levels:
    list_obj=list(f1["/%s" % (level)].keys())
    
    nstep=0
    for obj in list_obj:
        if "data" in obj:
            nstep = nstep + 1
    #
    list_steps.append(nstep)
#
print("steps = ", list_steps)
for i in range(len(list_steps)-1):
    lv = i+1
    if list_steps[i] < list_steps[i+1]:
        lvf1 = lv
        break
#
for i in range(len(list_steps)-1,0,-1):
    lv = i+1
    #print(i,lv,list_steps[i], list_steps[i-1])
    if list_steps[i] > list_steps[i-1]:
        lvf2 = lv
        break
#

print("lvf1,lvf2=",lvf1,lvf2)

n_lv = len(list_levels)
# lv_max = len(list_levels)

x=[]
y=[]
z=[]

list_data_lv = []
for level in list_levels:

    x.append(np.asarray(f1['/%s/x' % (level)]))
    y.append(np.asarray(f1['/%s/y' % (level)]))
    z.append(np.asarray(f1['/%s/z' % (level)]))

    list_obj=list(f1["/%s" % (level)].keys())

    list_obj.remove("x")
    list_obj.remove("y")
    list_obj.remove("z")
    # print(len(list_obj))
    list_data_lv.append(list_obj)

#list_obj=list(f1["/level%d" % (lv_max)].keys())
nstep_tot = len(list_data_lv[-1])
print(nstep_tot)

lv_max = n_lv

cmap_lv = cm.get_cmap('jet', 256)
norm_lv = colors.Normalize(vmin=1,vmax=lv_max)

for islice in range(nstep_tot):
    
    fig_xy = plt.figure(figsize=(size_fig, size_fig*0.95))
    ax_xy = fig_xy.add_subplot(111)
    ax_xy.set_xlabel("$x$ (cm)")
    ax_xy.set_ylabel("$y$ (cm)")
    ax_xy.set_xlim(-r_max,r_max)
    ax_xy.set_ylim(-r_max,r_max)
    ax_xy.set_aspect('equal')

    fig_yz = plt.figure(figsize=(size_fig, size_fig*0.95))
    ax_yz = fig_yz.add_subplot(111)
    ax_yz.set_xlabel("$y$ (cm)")
    ax_yz.set_ylabel("$z$ (cm)")
    ax_yz.set_xlim(-r_max,r_max)
    ax_yz.set_ylim(0,r_max)
    ax_yz.set_aspect('equal')

    fig_zx = plt.figure(figsize=(size_fig, size_fig*0.95))
    ax_zx = fig_zx.add_subplot(111)
    ax_zx.set_xlabel("$x$ (cm)")
    ax_zx.set_ylabel("$z$ (cm)")
    ax_zx.set_xlim(-r_max,r_max)
    ax_zx.set_ylim(0,r_max)
    ax_zx.set_aspect('equal')


    for ilv in range(n_lv):
        
        ax_xy.plot( [x[ilv][ 0],x[ilv][ 0]], [y[ilv][ 0],y[ilv][-1]], color="k", lw=1.0, ls="dotted" )
        ax_xy.plot( [x[ilv][-1],x[ilv][-1]], [y[ilv][ 0],y[ilv][-1]], color="k", lw=1.0, ls="dotted" )
        ax_xy.plot( [x[ilv][ 0],x[ilv][-1]], [y[ilv][ 0],y[ilv][ 0]], color="k", lw=1.0, ls="dotted" )
        ax_xy.plot( [x[ilv][ 0],x[ilv][-1]], [y[ilv][-1],y[ilv][-1]], color="k", lw=1.0, ls="dotted" )

        ax_yz.plot( [y[ilv][ 0],y[ilv][ 0]], [z[ilv][ 0],z[ilv][-1]], color="k", lw=1.0, ls="dotted" )
        ax_yz.plot( [y[ilv][-1],y[ilv][-1]], [z[ilv][ 0],z[ilv][-1]], color="k", lw=1.0, ls="dotted" )
        ax_yz.plot( [y[ilv][ 0],y[ilv][-1]], [z[ilv][ 0],z[ilv][ 0]], color="k", lw=1.0, ls="dotted" )
        ax_yz.plot( [y[ilv][ 0],y[ilv][-1]], [z[ilv][-1],z[ilv][-1]], color="k", lw=1.0, ls="dotted" )

        ax_zx.plot( [x[ilv][ 0],x[ilv][ 0]], [z[ilv][ 0],z[ilv][-1]], color="k", lw=1.0, ls="dotted" )
        ax_zx.plot( [x[ilv][-1],x[ilv][-1]], [z[ilv][ 0],z[ilv][-1]], color="k", lw=1.0, ls="dotted" )
        ax_zx.plot( [x[ilv][ 0],x[ilv][-1]], [z[ilv][ 0],z[ilv][ 0]], color="k", lw=1.0, ls="dotted" )
        ax_zx.plot( [x[ilv][ 0],x[ilv][-1]], [z[ilv][-1],z[ilv][-1]], color="k", lw=1.0, ls="dotted" )

    print("slice=",islice)
    text = "evolved in levels: "
    for ilv in range(n_lv):
        lv = ilv + 1

        if lv >= lvf2:
            nskip_lv = 1
        elif lvf1 <= lv:
            nskip_lv = 2**(lvf2-lv)
        else:
            nskip_lv = 2**(lvf2-lvf1)
        #
        if islice%nskip_lv==0:
            idat = int(islice/nskip_lv)
            group = list_data_lv[ilv][idat]
            group = "/level%d/" % (lv)+group

            t = np.asarray(f1[group+"/time"])[0]            
            Nevolved = np.asarray(f1[group+"/Nevolved"])[0]

            if Nevolved>0:
                text += "%d " % (lv)
                print(lv,idat,t,Nevolved)
                x_p = np.asarray(f1[group+"/x_evolved"])
                y_p = np.asarray(f1[group+"/y_evolved"])
                z_p = np.asarray(f1[group+"/z_evolved"])

                col = cmap_lv(norm_lv(lv))
                ax_xy.scatter(x_p, y_p, color=col, s=10, marker=".")
                ax_yz.scatter(y_p, z_p, color=col, s=10, marker=".")
                ax_zx.scatter(x_p, z_p, color=col, s=10, marker=".")
        #
    #print(text)
    ax_xy.text(1.0, 1.02, "$t=$%13.10f s" % (t), ha="right", va="top", transform=ax_xy.transAxes )
    ax_xy.text(1.0, 0.00, text, ha="left", va="top", transform=ax_xy.transAxes )
    ax_yz.text(1.0, 1.02, "$t=$%13.10f s" % (t), ha="right", va="top", transform=ax_yz.transAxes )
    ax_yz.text(1.0, 0.00, text, ha="left", va="top", transform=ax_xy.transAxes )
    ax_zx.text(1.0, 1.02, "$t=$%13.10f s" % (t), ha="right", va="top", transform=ax_zx.transAxes )
    ax_zx.text(1.0, 0.00, text, ha="left", va="top", transform=ax_xy.transAxes )
    #
    fn = dir_out + "/xy_%09d.png" % (islice)
    fig_xy.savefig(fn)
    plt.close(fig_xy)

    fn = dir_out + "/yz_%09d.png" % (islice)
    fig_yz.savefig(fn)
    plt.close(fig_yz)

    fn = dir_out + "/zx_%09d.png" % (islice)
    fig_zx.savefig(fn)
    plt.close(fig_zx)

    sys.exit()
    
