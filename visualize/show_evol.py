import numpy as np
# import h5py
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
from matplotlib.ticker import MaxNLocator,MultipleLocator,AutoMinorLocator
plt.switch_backend('agg')

# os.environ["PATH"] += os.pathsep + '/data/home/sfujibayashi/libs/texlive/2020/bin/x86_64-linux'

params={
    'font.size'           : 20.0     ,
    'font.family'         : 'DeJaVu Sans'  ,
#    'font.family'         : 'Times New Roman'  ,
    'xtick.major.size'    : 2        ,
    'xtick.major.width'   : 1.5      ,
    'xtick.labelsize'     : 20.0     ,
    'xtick.direction'     : 'in'     ,
    'ytick.major.size'    : 2        ,
    'ytick.major.width'   : 1.5      ,
    'ytick.labelsize'     : 20.0     ,
    'ytick.direction'     : 'in'     ,
    'xtick.major.pad'     : 2        ,
    'xtick.minor.pad'     : 2        ,
    'ytick.major.pad'     : 2        ,
    'ytick.minor.pad'     : 2        ,
    'axes.linewidth'      : 1.5      ,
    'text.usetex'         : False   }

plt.rcParams.update(params)


clight=2.99792458e10
r_uni = 4.81526840061178e4
rho_uni = 5.807834130139644e18
t_uni = 1.606200647186321e-6
msun = 1.989e33


ntraj_per_fig=1000
skip=100
alpha=0.4

#model = "DD2Tim326_Q4_M135_a75_0056_400m_B3e15_Hon_an20231126"; submodel = "ns9"

# model = "DD2Tim326_Q4_M135_a75_0056_400m_B3e15_Hon_an20231126"; submodel = "ns9"

# if len(submodel)>0:
#     submodel = "_"+submodel

# dir_read = "/raven/ptmp/shofu/"+model+"/Analysis_ptr/data"+submodel
# dir_out  = "/raven/ptmp/shofu/"+model+"/Analysis_ptr/fig"+submodel

parser = argparse.ArgumentParser()
parser.add_argument('--dir_read', required=False, default=".", type=str)
parser.add_argument('--dir_fig', required=False, default="./fig", type=str)

args = parser.parse_args()
dir_read = args.dir_read
dir_fig = args.dir_fig

#dir_read = "/sakura/ptmp/kiuchikn/DD2Tim326_135_135_0028_12.5mstg_B15.5_HLLD_CT_GS_reduced_lv14_to_lv13_run2/Analysis_ptr/data_3.0e8cm"
#dir_fig = "./fig_3.0e8cm"


fn="anim.sh"
with open(fn, mode="w") as f:
    f.write("ffmpeg -y -pattern_type glob -i '%s/rho_??????.png' -vf \"pad=ceil(iw/2)*2:ceil(ih/2)*2\" -pix_fmt yuv420p %s/rho.mp4\n" % (dir_fig, dir_fig))
    f.write("ffmpeg -y -pattern_type glob -i '%s/temp_??????.png' -vf \"pad=ceil(iw/2)*2:ceil(ih/2)*2\" -pix_fmt yuv420p %s/temp.mp4\n" % (dir_fig, dir_fig))
    f.write("ffmpeg -y -pattern_type glob -i '%s/entr_??????.png' -vf \"pad=ceil(iw/2)*2:ceil(ih/2)*2\" -pix_fmt yuv420p %s/entr.mp4\n" % (dir_fig, dir_fig))
    f.write("ffmpeg -y -pattern_type glob -i '%s/ye_??????.png' -vf \"pad=ceil(iw/2)*2:ceil(ih/2)*2\" -pix_fmt yuv420p %s/ye.mp4\n" % (dir_fig, dir_fig))
    f.write("ffmpeg -y -pattern_type glob -i '%s/r_??????.png' -vf \"pad=ceil(iw/2)*2:ceil(ih/2)*2\" -pix_fmt yuv420p %s/r.mp4\n" % (dir_fig, dir_fig))
    f.write("ffmpeg -y -pattern_type glob -i '%s/x_??????.png' -vf \"pad=ceil(iw/2)*2:ceil(ih/2)*2\" -pix_fmt yuv420p %s/x.mp4\n" % (dir_fig, dir_fig))
    f.write("ffmpeg -y -pattern_type glob -i '%s/y_??????.png' -vf \"pad=ceil(iw/2)*2:ceil(ih/2)*2\" -pix_fmt yuv420p %s/y.mp4\n" % (dir_fig, dir_fig))
    f.write("ffmpeg -y -pattern_type glob -i '%s/z_??????.png' -vf \"pad=ceil(iw/2)*2:ceil(ih/2)*2\" -pix_fmt yuv420p %s/z.mp4\n" % (dir_fig, dir_fig))
    f.write("ffmpeg -y -pattern_type glob -i '%s/vx_??????.png' -vf \"pad=ceil(iw/2)*2:ceil(ih/2)*2\" -pix_fmt yuv420p %s/vx.mp4\n" % (dir_fig, dir_fig))
    f.write("ffmpeg -y -pattern_type glob -i '%s/vy_??????.png' -vf \"pad=ceil(iw/2)*2:ceil(ih/2)*2\" -pix_fmt yuv420p %s/vy.mp4\n" % (dir_fig, dir_fig))
    f.write("ffmpeg -y -pattern_type glob -i '%s/vz_??????.png' -vf \"pad=ceil(iw/2)*2:ceil(ih/2)*2\" -pix_fmt yuv420p %s/vz.mp4\n" % (dir_fig, dir_fig))
    f.write("ffmpeg -y -pattern_type glob -i '%s/vr_??????.png' -vf \"pad=ceil(iw/2)*2:ceil(ih/2)*2\" -pix_fmt yuv420p %s/vr.mp4\n" % (dir_fig, dir_fig))
    f.write("ffmpeg -y -pattern_type glob -i '%s/rhoT_??????.png' -vf \"pad=ceil(iw/2)*2:ceil(ih/2)*2\" -pix_fmt yuv420p %s/rhoT.mp4\n" % (dir_fig, dir_fig))
    pass

if not os.path.exists(dir_fig):
    os.system('mkdir -p %s' % (dir_fig) )

list_data=sorted(glob.glob("%s/traj_????????.dat" % (dir_read)))
ntraj=len(list_data)

print(ntraj)
nfig = ntraj//ntraj_per_fig+1



for i in range(nfig):
    ip_start = 1+i*ntraj_per_fig
    ip_end   = np.min([1+(i+1)*ntraj_per_fig, ntraj+1])
    list_ip = range(ip_start,ip_end,skip)

    cmap_ip = plt.cm.get_cmap("jet")
    norm_ip = colors.Normalize(vmin=ip_start-1, vmax=ip_end)

    fig_rho = plt.figure(figsize=(10.0, 10.0*0.625))
    ax_rho  = fig_rho.add_subplot(111)
    #ax_rho.set_xlim(0.004,0.1)
    ax_rho.set_ylim(1e1,1e13)
    ax_rho.set_xscale("linear")
    ax_rho.set_yscale("log")
    ax_rho.set_xlabel("$t$ (s)")
    ax_rho.set_ylabel("$\\rho$ (g/cm$^3$)")
    ax_rho.xaxis.set_major_locator(MaxNLocator(nbins=4,min_n_ticks=2))
    ax_rho.xaxis.set_minor_locator(AutoMinorLocator(4))

    fig_temp = plt.figure(figsize=(10.0, 10.0*0.625))
    ax_temp  = fig_temp.add_subplot(111)
    #ax_temp.set_xlim(0.004,0.1)
    ax_temp.set_ylim(0.1,100.0)
    ax_temp.set_xscale("linear")
    ax_temp.set_yscale("log")
    ax_temp.set_xlabel("$t$ (s)")
    ax_temp.set_ylabel("$T$ (GK)")
    ax_temp.xaxis.set_major_locator(MaxNLocator(nbins=4,min_n_ticks=2))
    ax_temp.xaxis.set_minor_locator(AutoMinorLocator(4))

    fig_entr = plt.figure(figsize=(10.0, 10.0*0.625))
    ax_entr  = fig_entr.add_subplot(111)
    #ax_entr.set_xlim(0.004,0.1)
    ax_entr.set_ylim(1.0,300.0)
    ax_entr.set_xscale("linear")
    ax_entr.set_yscale("log")
    ax_entr.set_xlabel("$t$ (s)")
    ax_entr.set_ylabel("$s/k$ (nuc$^{-1}$)")
    ax_entr.xaxis.set_major_locator(MaxNLocator(nbins=4,min_n_ticks=2))
    ax_entr.xaxis.set_minor_locator(AutoMinorLocator(4))


    fig_ye = plt.figure(figsize=(10.0, 10.0*0.625))
    ax_ye  = fig_ye.add_subplot(111)
    #ax_ye.set_xlim(0.004,0.1)
    ax_ye.set_ylim(0.0,0.6)
    ax_ye.set_xscale("linear")
    ax_ye.set_yscale("linear")
    ax_ye.set_xlabel("$t$ (s)")
    ax_ye.set_ylabel("$Y_\\mathrm{e}$")
    ax_ye.xaxis.set_major_locator(MaxNLocator(nbins=4,min_n_ticks=2))
    ax_ye.xaxis.set_minor_locator(AutoMinorLocator(4))

    fig_r = plt.figure(figsize=(10.0, 10.0*0.625))
    ax_r  = fig_r.add_subplot(111)
    #ax_r.set_xlim(0.004,0.1)
    #ax_r.set_ylim(0.1,100.0)
    ax_r.set_xscale("linear")
    ax_r.set_yscale("log")
    ax_r.set_xlabel("$t$ (s)")
    ax_r.set_ylabel("$r$ (cm)")
    ax_r.xaxis.set_major_locator(MaxNLocator(nbins=4,min_n_ticks=2))
    ax_r.xaxis.set_minor_locator(AutoMinorLocator(4))

    fig_x = plt.figure(figsize=(10.0, 10.0*0.625))
    ax_x  = fig_x.add_subplot(111)
    #ax_x.set_xlim(0.004,0.1)
    #ax_x.set_ylim(0.1,100.0)
    ax_x.set_xscale("linear")
    ax_x.set_yscale("linear")
    ax_x.set_xlabel("$t$ (s)")
    ax_x.set_ylabel("$x$ (cm)")
    ax_x.xaxis.set_major_locator(MaxNLocator(nbins=4,min_n_ticks=2))
    ax_x.xaxis.set_minor_locator(AutoMinorLocator(4))

    fig_y = plt.figure(figsize=(10.0, 10.0*0.625))
    ax_y  = fig_y.add_subplot(111)
    #ax_y.set_xlim(0.004,0.1)
    #ax_y.set_ylim(0.1,100.0)
    ax_y.set_xscale("linear")
    ax_y.set_yscale("linear")
    ax_y.set_xlabel("$t$ (s)")
    ax_y.set_ylabel("$y$ (cm)")
    ax_y.xaxis.set_major_locator(MaxNLocator(nbins=4,min_n_ticks=2))
    ax_y.xaxis.set_minor_locator(AutoMinorLocator(4))

    fig_z = plt.figure(figsize=(10.0, 10.0*0.625))
    ax_z  = fig_z.add_subplot(111)
    #ax_z.set_xlim(0.004,0.1)
    #ax_z.set_ylim(0.1,100.0)
    ax_z.set_xscale("linear")
    ax_z.set_yscale("linear")
    ax_z.set_xlabel("$t$ (s)")
    ax_z.set_ylabel("$z$ (cm)")
    ax_z.xaxis.set_major_locator(MaxNLocator(nbins=4,min_n_ticks=2))
    ax_z.xaxis.set_minor_locator(AutoMinorLocator(4))

    fig_vx = plt.figure(figsize=(10.0, 10.0*0.625))
    ax_vx  = fig_vx.add_subplot(111)
    #ax_vx.set_xlim(0.004,0.1)
    #ax_vx.set_ylim(0.1,100.0)
    ax_vx.set_xscale("linear")
    ax_vx.set_yscale("linear")
    ax_vx.set_xlabel("$t$ (s)")
    ax_vx.set_ylabel("$v^x$ (cm/s)")
    ax_vx.xaxis.set_major_locator(MaxNLocator(nbins=4,min_n_ticks=2))
    ax_vx.xaxis.set_minor_locator(AutoMinorLocator(4))

    fig_vy = plt.figure(figsize=(10.0, 10.0*0.625))
    ax_vy  = fig_vy.add_subplot(111)
    #ax_vy.set_xlim(0.004,0.1)
    #ax_vy.set_ylim(0.1,100.0)
    ax_vy.set_xscale("linear")
    ax_vy.set_yscale("linear")
    ax_vy.set_xlabel("$t$ (s)")
    ax_vy.set_ylabel("$v^y$ (cm/s)")
    ax_vy.xaxis.set_major_locator(MaxNLocator(nbins=4,min_n_ticks=2))
    ax_vy.xaxis.set_minor_locator(AutoMinorLocator(4))

    fig_vz = plt.figure(figsize=(10.0, 10.0*0.625))
    ax_vz  = fig_vz.add_subplot(111)
    #ax_vz.set_xlim(0.004,0.1)
    #ax_vz.set_ylim(0.1,100.0)
    ax_vz.set_xscale("linear")
    ax_vz.set_yscale("linear")
    ax_vz.set_xlabel("$t$ (s)")
    ax_vz.set_ylabel("$v^z$ (cm/s)")
    ax_vz.xaxis.set_major_locator(MaxNLocator(nbins=4,min_n_ticks=2))
    ax_vz.xaxis.set_minor_locator(AutoMinorLocator(4))

    fig_vr = plt.figure(figsize=(10.0, 10.0*0.625))
    ax_vr  = fig_vr.add_subplot(111)
    #ax_vr.set_xlim(0.004,0.1)
    ax_vr.set_ylim(1e7,3e10)
    ax_vr.set_xscale("linear")
    ax_vr.set_yscale("log")
    ax_vr.set_xlabel("$t$ (s)")
    ax_vr.set_ylabel("$v^r$ (cm/s)")
    #ax_vr.xaxis.set_major_locator(MaxNLocator(nbins=4,min_n_ticks=2))
    #ax_vr.xaxis.set_minor_locator(AutoMinorLocator(4))

    fig_rhoT = plt.figure(figsize=(10.0, 10.0*0.625))
    ax_rhoT  = fig_rhoT.add_subplot(111)
    ax_rhoT.set_ylim(1e1,1e13)
    ax_rhoT.set_ylim(0.1,300.0)
    ax_rhoT.set_xscale("log")
    ax_rhoT.set_yscale("log")
    ax_rhoT.set_xlabel("$\\rho$ (g/cm$^3$)")
    ax_rhoT.set_ylabel("$T$ (GK)")

    for ip in list_ip:

        fn = dir_read + "/traj_%08d.dat" % (ip)
        
        print(ip)
        t  = np.loadtxt(fn,comments="#",usecols=(0))
        x  = np.loadtxt(fn,comments="#",usecols=(1))
        y  = np.loadtxt(fn,comments="#",usecols=(2))
        z  = np.loadtxt(fn,comments="#",usecols=(3))

        vx  = np.loadtxt(fn,comments="#",usecols=(4))
        vy  = np.loadtxt(fn,comments="#",usecols=(5))
        vz  = np.loadtxt(fn,comments="#",usecols=(6))

        rho = np.loadtxt(fn,comments="#",usecols=(7))
        temp= np.loadtxt(fn,comments="#",usecols=(8))
        ye = np.loadtxt(fn,comments="#",usecols=(9))
        entr = np.loadtxt(fn,comments="#",usecols=(10))
        
        r = np.sqrt(x**2+y**2+z**2)
        vr = (vx*x + vy*y + vz*z)/r

        col = cmap_ip(norm_ip(ip))

        ax_rho.plot(t,rho,color=col,alpha=alpha)
        ax_temp.plot(t,temp*1e-9,color=col,alpha=alpha)
        ax_entr.plot(t,entr,color=col,alpha=alpha)
        ax_ye.plot(t,ye,color=col,alpha=alpha)
        ax_r.plot(t,r,color=col,alpha=alpha)
        ax_x.plot(t,x,color=col,alpha=alpha)
        ax_y.plot(t,y,color=col,alpha=alpha)
        ax_z.plot(t,z,color=col,alpha=alpha)
        ax_vx.plot(t,vx,color=col,alpha=alpha)
        ax_vy.plot(t,vy,color=col,alpha=alpha)
        ax_vz.plot(t,vz,color=col,alpha=alpha)
        ax_vr.plot(t,vr,color=col,alpha=alpha)
        ax_rhoT.plot(rho,temp*1e-9,color=col,alpha=alpha)

        t_extend = np.linspace(t[-1], t[-1] + (t[-1]-t[0])*0.3, 10)
        rho_extend = rho[-1]*(t_extend/t[-1])**(-3)
        ax_rho.plot(t_extend,rho_extend,color=col,alpha=alpha, ls="dashed")

        pass


    s_map = cm.ScalarMappable(norm=norm_ip, cmap=cmap_ip)
    s_map.set_array([])
    cbar=fig_rho.colorbar(s_map,pad=0.01, ax=ax_rho); cbar.set_label("particle id",rotation = -90, labelpad= 30)
    cbar=fig_temp.colorbar(s_map,pad=0.01, ax=ax_temp); cbar.set_label("particle id",rotation = -90, labelpad= 30)
    cbar=fig_entr.colorbar(s_map,pad=0.01, ax=ax_entr); cbar.set_label("particle id",rotation = -90, labelpad= 30)
    cbar=fig_ye.colorbar(s_map,pad=0.01, ax=ax_ye); cbar.set_label("particle id",rotation = -90, labelpad= 30)
    cbar=fig_r.colorbar(s_map,pad=0.01, ax=ax_r); cbar.set_label("particle id",rotation = -90, labelpad= 30)
    cbar=fig_x.colorbar(s_map,pad=0.01, ax=ax_x); cbar.set_label("particle id",rotation = -90, labelpad= 30)
    cbar=fig_y.colorbar(s_map,pad=0.01, ax=ax_y); cbar.set_label("particle id",rotation = -90, labelpad= 30)
    cbar=fig_z.colorbar(s_map,pad=0.01, ax=ax_z); cbar.set_label("particle id",rotation = -90, labelpad= 30)
    cbar=fig_vx.colorbar(s_map,pad=0.01, ax=ax_vx); cbar.set_label("particle id",rotation = -90, labelpad= 30)
    cbar=fig_vy.colorbar(s_map,pad=0.01, ax=ax_vy); cbar.set_label("particle id",rotation = -90, labelpad= 30)
    cbar=fig_vz.colorbar(s_map,pad=0.01, ax=ax_vz); cbar.set_label("particle id",rotation = -90, labelpad= 30)
    cbar=fig_vr.colorbar(s_map,pad=0.01, ax=ax_vr); cbar.set_label("particle id",rotation = -90, labelpad= 30)
    cbar=fig_rhoT.colorbar(s_map,pad=0.01, ax=ax_rhoT); cbar.set_label("particle id",rotation = -90, labelpad= 30)
    #cbar.set_ticks([6.0/6.0, 7.0/6.0, 8.0/6.0, 9.0/6.0, 10.0/6.0])
    #cbar.ax.set_yticklabels(["1", "7/6", "4/3", "3/2", "5/3"])

    fig_rho.subplots_adjust(left=0.15,bottom=0.15,right=0.95, top=0.95)
    fig_temp.subplots_adjust(left=0.15,bottom=0.15,right=0.95, top=0.95)
    fig_entr.subplots_adjust(left=0.15,bottom=0.15,right=0.95, top=0.95)
    fig_ye.subplots_adjust(left=0.15,bottom=0.15,right=0.95, top=0.95)
    fig_r.subplots_adjust(left=0.15,bottom=0.15,right=0.95, top=0.95)
    fig_x.subplots_adjust(left=0.15,bottom=0.15,right=0.95, top=0.95)
    fig_y.subplots_adjust(left=0.15,bottom=0.15,right=0.95, top=0.95)
    fig_z.subplots_adjust(left=0.15,bottom=0.15,right=0.95, top=0.95)
    fig_vx.subplots_adjust(left=0.15,bottom=0.15,right=0.95, top=0.95)
    fig_vy.subplots_adjust(left=0.15,bottom=0.15,right=0.95, top=0.95)
    fig_vz.subplots_adjust(left=0.15,bottom=0.15,right=0.95, top=0.95)
    fig_vr.subplots_adjust(left=0.15,bottom=0.15,right=0.95, top=0.95)
    fig_rhoT.subplots_adjust(left=0.15,bottom=0.15,right=0.95, top=0.95)

    fig_rho.savefig(dir_fig + "/rho_%06d.png" % (i))
    fig_temp.savefig(dir_fig + "/temp_%06d.png" % (i))
    fig_entr.savefig(dir_fig + "/entr_%06d.png" % (i))
    fig_ye.savefig(dir_fig + "/ye_%06d.png" % (i))
    fig_r.savefig(dir_fig + "/r_%06d.png" % (i))
    fig_x.savefig(dir_fig + "/x_%06d.png" % (i))
    fig_y.savefig(dir_fig + "/y_%06d.png" % (i))
    fig_z.savefig(dir_fig + "/z_%06d.png" % (i))
    fig_vx.savefig(dir_fig + "/vx_%06d.png" % (i))
    fig_vy.savefig(dir_fig + "/vy_%06d.png" % (i))
    fig_vz.savefig(dir_fig + "/vz_%06d.png" % (i))
    fig_vr.savefig(dir_fig + "/vr_%06d.png" % (i))
    fig_rhoT.savefig(dir_fig + "/rhoT_%06d.png" % (i))

    plt.close(fig_rho)
    plt.close(fig_temp)
    plt.close(fig_entr)
    plt.close(fig_ye)
    plt.close(fig_r)
    plt.close(fig_x)
    plt.close(fig_y)
    plt.close(fig_z)
    plt.close(fig_vx)
    plt.close(fig_vy)
    plt.close(fig_vz)
    plt.close(fig_vr)
    plt.close(fig_rhoT)

    pass
