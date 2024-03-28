import numpy as np
# import h5py
import sys
#import subprocess
import os
import glob
from scipy.interpolate import interp2d

import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as colors
import matplotlib.ticker as ticker
plt.switch_backend('agg')
os.environ["PATH"] += os.pathsep + '/data/home/sfujibayashi/libs/texlive/2020/bin/x86_64-linux'

params={
    'font.size'           : 24.0     ,
    'font.family'         : 'DeJaVu Sans'  ,
#    'font.family'         : 'Times New Roman'  ,
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
    'axes.linewidth'      : 1.5      ,
    'text.usetex'         : True    }

plt.rcParams.update(params)


clight=2.99792458e10
r_uni = 4.81526840061178e4
rho_uni = 5.807834130139644e18
t_uni = 1.606200647186321e-6
msun = 1.989e33

dir_base="/scratch/sfujibayashi/Particle_trace_data/data_3D_BHNSDD2q4_hr/small_set"

dir_read = dir_base + "/data_dyn_unified"
dir_out  = dir_base + "/fig"
label="DD2q4_high"

t_merge=9.35e-3
os.system('mkdir -p %s' % (dir_out) )


data_traj=np.loadtxt(dir_read + "/ana_traj.dat" ,comments='#')


#print(data_traj[:,0])
id_tracer = data_traj[:,0].astype(int)
condition_ut1 = data_traj[:,15]
condition_hut = data_traj[:,16]
number_tracer = len(id_tracer)
time_boundary = data_traj[:,3]
mass_tracer = data_traj[:,2]

mass_tracer_total = np.sum(mass_tracer)

fig1 = plt.figure(figsize=(10.0, 10.0*0.625))
ax1  = fig1.add_subplot(111)
ax1.hist(condition_ut1, weights=mass_tracer/mass_tracer_total,range=([-0.07,0.03]),bins=30,alpha=0.5,label="$u_t + 1$")
ax1.hist(condition_hut, weights=mass_tracer/mass_tracer_total,range=([-0.07,0.03]),bins=30,alpha=0.5,label="$hu_t + h_\mathrm{min}$")

fig1.subplots_adjust(left=0.15,bottom=0.15,right=0.95, top=0.95)
ax1.legend(fontsize=16,frameon=False)

ax1.set_xlabel("$hu_t+h_\mathrm{min}$ or $u_1 +1$")
ax1.set_ylabel("$\Delta M/M_\\mathrm{tot}$ (cm)")

fig1.savefig("hist_ut1.png")


list_data=sorted(glob.glob("%s/traj_????????.dat" % (dir_read)))
number_tracer=len(list_data)

# fig1 = plt.figure(figsize=(10.0, 10.0*0.625))
# ax1  = fig1.add_subplot(111)
# fig2 = plt.figure(figsize=(10.0, 10.0*0.625))
# ax2  = fig2.add_subplot(111)
# fig3 = plt.figure(figsize=(10.0, 10.0*0.625))
# ax3  = fig3.add_subplot(111)

fig_vel_t = plt.figure(figsize=(10.0, 10.0*0.625))
ax_vel_t  = fig_vel_t.add_subplot(111)

fig_vel_r = plt.figure(figsize=(10.0, 10.0*0.625))
ax_vel_r  = fig_vel_r.add_subplot(111)

mass_target=0.0
n_haffected=0
for i in range(number_tracer):
    ip = id_tracer[i]
    fn=dir_read + "/traj_%08i.dat" % (ip)
    data = np.loadtxt(fn,comments='#')
    t_tracer = data[:,0]
    x_tracer = data[:,1]
    y_tracer = data[:,2]
    z_tracer = data[:,3]
    vx_tracer = data[:,4]
    vy_tracer = data[:,5]
    vz_tracer = data[:,6]
    
    rho_tracer=data[:,7]
    tem_tracer=data[:,8]
    ye_tracer =data[:,9]
    sen_tracer =data[:,10]

    ut1_tracer = data[:,13]
    hut_tracer = data[:,14]
    
    r_tracer = np.sqrt(x_tracer**2+y_tracer**2+z_tracer**2)
    
    vr_tracer = (vx_tracer*x_tracer + vy_tracer*y_tracer + vz_tracer*z_tracer)/r_tracer
    
    #it_diag = np.argmin( np.abs(t_tracer-t_merge - 10e-3) )
    it_diag=len(t_tracer)-1
    
    d_condition = ut1_tracer[it_diag] - hut_tracer[it_diag]
    
    # alpha = max([0.1,min([1.0,mass_tracer[i]/2e31 * 10.0])])
    # if r_tracer[0] < 5e7 and d_condition > 0.007:
    #     n_haffected = n_haffected + 1
    #     if ut1_tracer[it_diag] >0.0 and hut_tracer[it_diag] <0.0:
    #         mass_target = mass_target + mass_tracer[i]
    #         ax1.plot(t_tracer-t_merge, r_tracer ,color="k", ls = "solid", alpha = alpha, linewidth = 1.0)
    #         ax2.plot(t_tracer-t_merge, hut_tracer ,color="k", ls = "solid", alpha = alpha, linewidth = 1.0)
    #         ax2.plot(t_tracer-t_merge, ut1_tracer ,color="k", ls = "dashed", alpha = alpha, linewidth = 1.0)
    #         print( ("%6i"+"%12.4e"*5) % (ip,r_tracer[it_diag],vr_tracer[it_diag],ut1_tracer[it_diag],hut_tracer[it_diag], mass_tracer[i]) )
    #     elif ut1_tracer[it_diag] < 0.0 and hut_tracer[it_diag] <0.0:
    #         ax1.plot(t_tracer-t_merge, r_tracer ,color="b", ls = "solid", alpha = alpha, linewidth = 1.0)
    #         ax3.plot(t_tracer-t_merge, hut_tracer ,color="b", ls = "solid", alpha = alpha, linewidth = 1.0)
    #         ax3.plot(t_tracer-t_merge, ut1_tracer ,color="b", ls = "dashed", alpha = alpha, linewidth = 1.0)
        
        #ax2.plot(t_tracer-t_merge, hut_tracer ,color="b", ls = "dashed", alpha = 0.1, linewidth = 1.0)
    if ye_tracer[-1]<0.08:


        if np.amax(rho_tracer)>4e11:

            if condition_ut1[i]<0.0:
                ax_vel_r.plot(r_tracer, vr_tracer/clight ,color="k", ls = "solid", alpha = 0.2, linewidth = 1.0)
                ax_vel_t.plot(t_tracer-t_merge, vr_tracer/clight ,color="k", ls = "solid", alpha = 0.2, linewidth = 1.0)
            else:
                ax_vel_r.plot(r_tracer, vr_tracer/clight ,color="b", ls = "solid", alpha = 0.2, linewidth = 1.0)
                ax_vel_t.plot(t_tracer-t_merge, vr_tracer/clight ,color="b", ls = "solid", alpha = 0.2, linewidth = 1.0)

            for it in range(len(rho_tracer)-1):
                if rho_tracer[it] >= 4e11 and 4e11 > rho_tracer[it+1]:
                    it_drip = it
                    break
            # for it in range(it_drip-10,it_drip+10):
            #     print( ("%5i"+"%12.4e"*3) % (it,t_tracer[it],rho_tracer[it],vr_tracer[it]/clight))
            t1 = (4e11-rho_tracer[it_drip])/(rho_tracer[it_drip+1]-rho_tracer[it_drip])
            #t1 = (np.log10(4e11)-np.log10(rho_tracer[it_drip]))/(np.log10(rho_tracer[it_drip+1])-np.log10(rho_tracer[it_drip]))
            t0 = 1.0-t1
            
            t_drip = t0* t_tracer[it_drip]+t1* t_tracer[it_drip+1]
            r_drip = t0* r_tracer[it_drip]+t1* r_tracer[it_drip+1]
            vr_drip= t0*vr_tracer[it_drip]+t1*vr_tracer[it_drip+1]
            
            #if vr_drip/clight<0.1:
                #ax_vel_t.plot(t_tracer-t_merge, vr_tracer/clight ,color="r", ls = "solid", alpha = 0.2, linewidth = 1.0)
            if condition_ut1[i]<0.0:
                ax_vel_t.scatter(t_drip-t_merge, vr_drip/clight, color="k",s=10)
                ax_vel_r.scatter(r_drip, vr_drip/clight, color="k",s=10)
            else:
                #if vr_drip/clight>0.2:
                ax_vel_t.scatter(t_drip-t_merge, vr_drip/clight, color="b",s=10)
                ax_vel_r.scatter(r_drip, vr_drip/clight, color="b",s=10)
            # print(ip,it_drip)
            # print(("%12.4e"*3) % (t1,rho_tracer[it_drip],rho_tracer[it_drip+1]))
            # print(("%12.4e"*3) % (vr_tracer[it_drip]/clight,vr_drip/clight,vr_tracer[it_drip+1]/clight))
            
            # break
        
x_tmp=[-10,-9]
ax_vel_r.plot(x_tmp,x_tmp ,color="b", ls = "solid", alpha = 1, linewidth = 1.0,label="$u_t>-1$ (still bound)")
ax_vel_t.plot(x_tmp,x_tmp ,color="b", ls = "solid", alpha = 1, linewidth = 1.0,label="$u_t>-1$ (still bound)")

# print(mass_target,mass_tracer_total,n_haffected)
td=-3e-3;tu=7e-3;scalex="linear"

# fig1.subplots_adjust(left=0.15,bottom=0.15,right=0.95, top=0.95)
# ax1.legend(fontsize=16,frameon=False)
# ax1.set_xlim(td,tu)
# ax1.set_ylim(1e6,1e9)
# ax1.set_xscale(scalex)
# ax1.set_yscale('log')
# ax1.tick_params(axis='both',width=1.5,length=6,which='major',top=True,right=True,pad=10.0)
# ax1.tick_params(axis='both',width=1.5,length=3,which='minor',top=True,right=True)
# ax1.set_xlabel("$t-t_\\mathrm{merge}$ (s)")
# ax1.set_ylabel("$r$ (cm)")
# #ax1.xaxis.set_major_locator(ticker.MultipleLocator(1.0))
# #ax1.xaxis.set_minor_locator(ticker.MultipleLocator(0.2))
# fig1.savefig("tr_x%s.png" % (scalex))

# fig2.subplots_adjust(left=0.15,bottom=0.15,right=0.95, top=0.95)
# ax2.legend(fontsize=16,frameon=False)
# ax2.set_xlim(td,tu)
# ax2.set_ylim(-0.1,+0.1)
# ax2.set_xscale(scalex)
# #ax2.set_yscale('log')
# ax2.tick_params(axis='both',width=1.5,length=6,which='major',top=True,right=True,pad=10.0)
# ax2.tick_params(axis='both',width=1.5,length=3,which='minor',top=True,right=True)
# ax2.set_xlabel("$t-t_\\mathrm{merge}$ (s)")
# ax2.set_ylabel("$h u_t + h_\mathrm{min}$ or $u_t + 1$")
# #ax2.xaxis.set_major_locator(ticker.MultipleLocator(1.0))
# #ax2.xaxis.set_minor_locator(ticker.MultipleLocator(0.2))
# fig2.savefig("tr_hut_x%s.png" % (scalex))

# fig3.subplots_adjust(left=0.15,bottom=0.15,right=0.95, top=0.95)
# ax3.legend(fontsize=16,frameon=False)
# ax3.set_xlim(td,tu)
# ax3.set_ylim(-0.1,+0.1)
# ax3.set_xscale(scalex)
# #ax3.set_yscale('log')
# ax3.tick_params(axis='both',width=1.5,length=6,which='major',top=True,right=True,pad=10.0)
# ax3.tick_params(axis='both',width=1.5,length=3,which='minor',top=True,right=True)
# ax3.set_xlabel("$t-t_\\mathrm{merge}$ (s)")
# ax3.set_ylabel("$h u_t + h_\mathrm{min}$ or $u_t + 1$")
# #ax3.xaxis.set_major_locator(ticker.MultipleLocator(1.0))
# #ax3.xaxis.set_minor_locator(ticker.MultipleLocator(0.2))
# fig3.savefig("tr_hut_2_x%s.png" % (scalex))

# td=0.0001;tu=15e-3;scalex="log"
# fig1.subplots_adjust(left=0.15,bottom=0.15,right=0.95, top=0.95)
# ax1.legend(fontsize=16,frameon=False)
# ax1.set_xlim(td,tu)
# ax1.set_ylim(1e6,1e9)
# ax1.set_xscale(scalex)
# ax1.set_yscale('log')
# ax1.tick_params(axis='both',width=1.5,length=6,which='major',top=True,right=True,pad=10.0)
# ax1.tick_params(axis='both',width=1.5,length=3,which='minor',top=True,right=True)
# ax1.set_xlabel("$t-t_\\mathrm{merge}$ (s)")
# ax1.set_ylabel("$r$ (cm)")
# #ax1.xaxis.set_major_locator(ticker.MultipleLocator(1.0))
# #ax1.xaxis.set_minor_locator(ticker.MultipleLocator(0.2))
# fig1.savefig("tr_x%s.png" % (scalex))


fig_vel_t.subplots_adjust(left=0.15,bottom=0.15,right=0.95, top=0.95)
ax_vel_t.legend(fontsize=16,frameon=False)
ax_vel_t.set_xlim(td,tu)
ax_vel_t.set_ylim(-0.1,0.4)
ax_vel_t.set_xscale(scalex)
#ax_vel_t.set_yscale('log')
ax_vel_t.tick_params(axis='both',width=1.5,length=6,which='major',top=True,right=True,pad=10.0)
ax_vel_t.tick_params(axis='both',width=1.5,length=3,which='minor',top=True,right=True)
ax_vel_t.set_xlabel("$t-t_\\mathrm{merge}$ (s)")
ax_vel_t.set_ylabel("$v^r/c$")
#ax_vel_t.xaxis.set_major_locator(ticker.MultipleLocator(1.0))
#ax_vel_t.xaxis.set_minor_locator(ticker.MultipleLocator(0.2))
fig_vel_t.savefig("tr_vel_t_x%s_%s.png" % (scalex,label))

fig_vel_r.subplots_adjust(left=0.15,bottom=0.15,right=0.95, top=0.95)
ax_vel_r.legend(fontsize=16,frameon=False)
ax_vel_r.set_xlim(1e7,1e9)
ax_vel_r.set_ylim(-0.1,0.4)
ax_vel_r.set_xscale("log")
#ax_vel_r.set_yscale('log')
ax_vel_r.tick_params(axis='both',width=1.5,length=6,which='major',top=True,right=True,pad=10.0)
ax_vel_r.tick_params(axis='both',width=1.5,length=3,which='minor',top=True,right=True)
ax_vel_r.set_xlabel("$r$ (cm)")
ax_vel_r.set_ylabel("$v^r/c$")
#ax_vel_r.xaxis.set_major_locator(ticker.MultipleLocator(1.0))
#ax_vel_r.xaxis.set_minor_locator(ticker.MultipleLocator(0.2))
fig_vel_r.savefig("tr_vel_r_xlog_%s.png" % (label))

