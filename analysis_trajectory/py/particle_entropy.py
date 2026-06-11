import numpy as np
import h5py
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

### EOS
fn = "/scratch/sfujibayashi/EOS/EOS_Hempel_DD2Tim_TF_326.h5"
f0 = h5py.File(fn, 'r')
tem_e = np.asarray(f0['/tem'])
rho_e = np.asarray(f0['/rho'])
ye_e  = np.asarray(f0['/ye'])
aa_e  = np.asarray(f0['/aa'])
zz_e  = np.asarray(f0['/zz'])
xA_e  = np.asarray(f0['/xA'])

# take logalism
rho_e = np.log10(rho_e)
tem_e = np.log10(tem_e)
#eps_e = np.log10(eps_e/clight**2+1.0)
#pres_e = np.log10(pres_e/clight**2/rho_uni)
#cs_e  = np.log10(cs_e)

ied=0
ieu=len(tem_e)-1
ked=0
keu=len(rho_e)-1
jed=0
jeu=len(ye_e)-1

tem_e_min = tem_e[ied]
rho_e_min = rho_e[ked]
ye_e_min = ye_e[jed]
dtem = tem_e[1]-tem_e[0]
drho = rho_e[1]-rho_e[0]
dye = ye_e[1]-ye_e[0]
dtemi=1.0/dtem
drhoi=1.0/drho
dyei=1.0/dye
###

clight=2.99792458e10
r_uni = 4.81526840061178e4
rho_uni = 5.807834130139644e18
t_uni = 1.606200647186321e-6
msun = 1.989e33
k2mev = 1.0/1.160445e10

#dir_read = "/scratch/sfujibayashi/Particle_trace_data/data_3D_BHNSDD2q6_10ms/data"
#dir_out  = "/scratch/sfujibayashi/Particle_trace_data/data_3D_BHNSDD2q6_10ms/fig"
dir_read = "/scratch/sfujibayashi/Particle_trace_data/data_3D_BHNSDD2q6/data_post"
dir_out  = "/scratch/sfujibayashi/Particle_trace_data/data_3D_BHNSDD2q6/fig_post"

t_merge=9.35e-3
os.system('mkdir -p %s' % (dir_out) )


data_traj=np.loadtxt(dir_read + "/ana_traj.dat" ,comments='#')

#print(data_traj[:,0])
id_tracer = np.around(data_traj[:,0]).astype(int)
condition_ut1 = data_traj[:,15]
condition_hut = data_traj[:,16]
number_tracer = len(id_tracer)
time_final_tracer = data_traj[:,3]
mass_tracer = data_traj[:,2]
ye_final_tracer = data_traj[:,9]
temp_max_tracer =  data_traj[:,14]
mass_tracer_total = np.sum(mass_tracer)

#print(np.amax(temp_max_tracer)*k2mev)
#sys.exit()

# fig1 = plt.figure(figsize=(10.0, 10.0*0.625))
# ax1  = fig1.add_subplot(111)
# ax1.hist(condition_ut1, weights=mass_tracer/mass_tracer_total,range=([-0.07,0.03]),bins=30,alpha=0.5,label="$u_t + 1$")
# ax1.hist(condition_hut, weights=mass_tracer/mass_tracer_total,range=([-0.07,0.03]),bins=30,alpha=0.5,label="$hu_t + h_\mathrm{min}$")

# fig1.subplots_adjust(left=0.15,bottom=0.15,right=0.95, top=0.95)
# ax1.legend(fontsize=16,frameon=False)

# ax1.set_xlabel("$hu_t+h_\mathrm{min}$ or $u_1 +1$")
# ax1.set_ylabel("$\Delta M/M_\\mathrm{tot}$ (cm)")

# fig1.savefig("hist_ut1.png")


list_data=sorted(glob.glob("%s/traj_????????.dat" % (dir_read)))
number_tracer=len(list_data)

fig1 = plt.figure(figsize=(10.0, 10.0*0.625))
ax1  = fig1.add_subplot(111)

n_tracer=0
for i in range(0,number_tracer,5):

    ip = id_tracer[i]

    print(ip,time_final_tracer[i], temp_max_tracer[i]*k2mev)
    if time_final_tracer[i]>0.5 and temp_max_tracer[i]*k2mev>5.0:
        n_tracer = n_tracer + 1

        print(ip)
        
        fn=dir_read + "/traj_%08i.dat" % (ip)
        data = np.loadtxt(fn,comments='#')
        t_tracer = data[:,0]
        # x_tracer = data[:,1]
        # y_tracer = data[:,2]
        # z_tracer = data[:,3]
        # vx_tracer = data[:,4]
        # vy_tracer = data[:,5]
        # vz_tracer = data[:,6]
        
        # rho_tracer=data[:,7]
        tem_tracer=data[:,8]
        ye_tracer =data[:,9]
        # sen_tracer =data[:,10]

        it_plot = 0
        it_plot = np.argmin(np.abs(t_tracer-0.03))
        # for it in range(len(t_tracer)-1,2,-1):
        #     # print((tem_tracer[it]-tem_tracer[it-1])/tem_tracer[it])
        #     if np.abs((tem_tracer[it]-tem_tracer[it-1])/tem_tracer[it]) > 0.5:
        #         it_plot = it
        #         break
        # sys.exit()
        color="k"
        zorder=0
        alpha=0.3
        ax1.plot(tem_tracer[it_plot:]*k2mev, ye_tracer[it_plot:]   ,color=color, ls = "solid", alpha = alpha, linewidth = 1.0, zorder=zorder)

print(n_tracer)    
#print(mass_tracer_total,mass_ye01,mass_ye01_sjump,n_tracer,mass_ye01_yjump)
tempd=5;tempu=0.1;scalex="log"

fig1.subplots_adjust(left=0.15,bottom=0.15,right=0.95, top=0.95)
ax1.legend(fontsize=16,frameon=False)
ax1.set_xlim(tempd,tempu)
ax1.set_ylim(0.0,0.5)
ax1.set_xscale(scalex)
ax1.set_yscale('linear')
ax1.tick_params(axis='both',width=1.5,length=6,which='major',top=True,right=True,pad=10.0)
ax1.tick_params(axis='both',width=1.5,length=3,which='minor',top=True,right=True)
ax1.set_xlabel("$k_\\mathrm{B}T$ (MeV)")
ax1.set_ylabel("$Y_\\mathrm{e}$")
#ax1.xaxis.set_major_locator(ticker.MultipleLocator(1.0))
#ax1.xaxis.set_minor_locator(ticker.MultipleLocator(0.2))
fig1.savefig("temp_Ye_Q6.png")

