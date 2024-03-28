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

dir_read = "/scratch/sfujibayashi/Particle_trace_data/data_3D_BHNSDD2q6_10ms/data"
dir_out  = "/scratch/sfujibayashi/Particle_trace_data/data_3D_BHNSDD2q6_10ms/fig"

t_merge=9.35e-3
os.system('mkdir -p %s' % (dir_out) )


print(dir_read + "/ana_traj.dat")

data_traj=np.loadtxt(dir_read + "/ana_traj.dat" ,comments='#')

#print(data_traj[:,0])
id_tracer = np.around(data_traj[:,0])
condition_ut = data_traj[:,15]
condition_hut = data_traj[:,16]
number_tracer = len(id_tracer)
time_boundary = data_traj[:,3]
mass_tracer = data_traj[:,2]

mass_tracer_total = np.sum(mass_tracer)
mass_hut_but_ut = 0.0

tracer_hut_but_ut = []
tracer_hut_and_ut = []


for i in range(number_tracer):
    if condition_hut[i] < 0.0 and condition_ut[i] > 0.0:
        # print( ("%6i" + "%12.4e"*4) % (id_tracer[i],time_boundary[i],condition_ut[i],condition_hut[i],mass_tracer[i]) )
        tracer_hut_but_ut.append(i)

    if condition_hut[i] < 0.0 and condition_ut[i] < 0.0:
        tracer_hut_and_ut.append(i)


    list_data=sorted(glob.glob("%s/traj_????????.dat" % (dir_read)))
    number_tracer=len(list_data)
    # sys.exit()
    

        # mass_hut_but_ut = mass_hut_but_ut + mass_tracer[i]

# print(mass_hut_but_ut/msun,mass_tracer_total/msun)

fig1 = plt.figure(figsize=(10.0, 10.0*0.625))
ax1  = fig1.add_subplot(111)
# fig2 = plt.figure(figsize=(10.0, 10.0*0.625))
# ax2  = fig2.add_subplot(111)

for i in tracer_hut_and_ut[::10]:
    ip = id_tracer[i]
    if time_boundary[i] < 0.5:
        print(ip)
        fn=dir_read + "/traj_%08i.dat" % (ip)
        data = np.loadtxt(fn,comments='#')
        t_tracer = data[:,0]
        x_tracer = data[:,1]
        y_tracer = data[:,2]
        z_tracer = data[:,3]
        
        rho_tracer=data[:,7]
        tem_tracer=data[:,8]
        ye_tracer =data[:,9]
        
        r_tracer = np.sqrt(x_tracer**2+y_tracer**2+z_tracer**2)
        ax1.plot(t_tracer-t_merge, r_tracer ,color="r", ls = "solid", alpha = 0.2, linewidth = 1.0)
        #ax2.plot(t_tracer-t_merge, r_tracer ,color="r", ls = "solid", alpha = 0.2, linewidth = 1.0)

mass_target = 0.0
for i in tracer_hut_but_ut:
    ip = id_tracer[i]
    if time_boundary[i] < 0.5 or 1==1:
        fn=dir_read + "/traj_%08i.dat" % (ip)
        data = np.loadtxt(fn,comments='#')
        t_tracer = data[:,0]
        x_tracer = data[:,1]
        y_tracer = data[:,2]
        z_tracer = data[:,3]
        vx_tracer = data[:,1]
        vy_tracer = data[:,2]
        vz_tracer = data[:,3]
        
        rho_tracer=data[:,7]
        tem_tracer=data[:,8]
        ye_tracer =data[:,9]
        
        r_tracer = np.sqrt(x_tracer**2+y_tracer**2+z_tracer**2)

        vr_tracer = (vx_tracer*x_tracer + vy_tracer*y_tracer + vz_tracer*z_tracer)/r_tracer

        it_10pm = np.argmin( np.abs(t_tracer-t_merge - 10e-3) )
        print(ip,r_tracer[it_10pm],vr_tracer[it_10pm])
        if r_tracer[it_10pm] > 2e7:
            mass_target = mass_target + mass_tracer[i]
            ax1.plot(t_tracer-t_merge, r_tracer ,color="k", ls = "solid", alpha = 0.1, linewidth = 1.0)
        else:
            ax1.plot(t_tracer-t_merge, r_tracer ,color="b", ls = "solid", alpha = 0.1, linewidth = 1.0)
        #ax2.plot(t_tracer-t_merge, r_tracer ,color="k", ls = "solid", alpha = 0.2, linewidth = 1.0)

print(mass_target/msun)
        
td=-t_merge;tu=2e-2;scalex="linear"

fig1.subplots_adjust(left=0.15,bottom=0.15,right=0.95, top=0.95)
ax1.legend(fontsize=16,frameon=False)
ax1.set_xlim(td,tu)
ax1.set_ylim(1e6,1e9)
ax1.set_xscale(scalex)
ax1.set_yscale('log')
ax1.tick_params(axis='both',width=1.5,length=6,which='major',top=True,right=True,pad=10.0)
ax1.tick_params(axis='both',width=1.5,length=3,which='minor',top=True,right=True)
ax1.set_xlabel("$t-t_\\mathrm{merge}$ (s)")
ax1.set_ylabel("$r$ (cm)")
#ax1.xaxis.set_major_locator(ticker.MultipleLocator(1.0))
#ax1.xaxis.set_minor_locator(ticker.MultipleLocator(0.2))
fig1.savefig("tr_x%s.png" % (scalex))

td=0.0001;tu=1.4;scalex="log"
fig1.subplots_adjust(left=0.15,bottom=0.15,right=0.95, top=0.95)
ax1.legend(fontsize=16,frameon=False)
ax1.set_xlim(td,tu)
ax1.set_ylim(1e6,1e9)
ax1.set_xscale(scalex)
ax1.set_yscale('log')
ax1.tick_params(axis='both',width=1.5,length=6,which='major',top=True,right=True,pad=10.0)
ax1.tick_params(axis='both',width=1.5,length=3,which='minor',top=True,right=True)
ax1.set_xlabel("$t-t_\\mathrm{merge}$ (s)")
ax1.set_ylabel("$r$ (cm)")
#ax1.xaxis.set_major_locator(ticker.MultipleLocator(1.0))
#ax1.xaxis.set_minor_locator(ticker.MultipleLocator(0.2))
fig1.savefig("tr_x%s.png" % (scalex))
