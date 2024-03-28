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

#dir_read = "/scratch/sfujibayashi/Particle_trace_data/data_3D_BHNSDD2q6_10ms/data"
#dir_out  = "/scratch/sfujibayashi/Particle_trace_data/data_3D_BHNSDD2q6_10ms/fig"
dir_read = "/scratch/sfujibayashi/Particle_trace_data/data_3D_BHNSDD2q6/data"
dir_out  = "/scratch/sfujibayashi/Particle_trace_data/data_3D_BHNSDD2q6/fig"

t_merge=9.35e-3
os.system('mkdir -p %s' % (dir_out) )


data_traj=np.loadtxt(dir_read + "/ana_traj.dat" ,comments='#')

#print(data_traj[:,0])
id_tracer = np.around(data_traj[:,0])
condition_ut1 = data_traj[:,15]
condition_hut = data_traj[:,16]
number_tracer = len(id_tracer)
time_boundary = data_traj[:,3]
mass_tracer = data_traj[:,2]
ye_last_tracer = data_traj[:,9]
mass_tracer_total = np.sum(mass_tracer)

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
fig2 = plt.figure(figsize=(10.0, 10.0*0.625))
ax2  = fig2.add_subplot(111)
fig3 = plt.figure(figsize=(10.0, 10.0*0.625))
ax3  = fig3.add_subplot(111)
fig4 = plt.figure(figsize=(10.0, 10.0*0.625))
ax4  = fig4.add_subplot(111)
fig5 = plt.figure(figsize=(10.0, 10.0*0.625))
ax5  = fig5.add_subplot(111)
fig6 = plt.figure(figsize=(10.0, 10.0*0.625))
ax6  = fig6.add_subplot(111)
fig7 = plt.figure(figsize=(10.0, 10.0*0.625))
ax7  = fig7.add_subplot(111)

tmax_ana = 1.0e10

mass_ye01 = 0.0
mass_ye01_sjump=0.0
mass_ye01_yjump=0.0

n_tracer=0
for i in range(number_tracer):

    ip = id_tracer[i]
    dsen_max = 0.0
    if ye_last_tracer[i]<0.1:

        mass_ye01 = mass_ye01 + mass_tracer[i]
        
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
        # sen_tracer =np.maximum(data[:,10],1e-2)
        sen_tracer =data[:,10]
        
        #ut1_tracer = data[:,13]
        #hut_tracer = data[:,14]
        
        r_tracer = np.sqrt(x_tracer**2+y_tracer**2+z_tracer**2)
        
        vr_tracer = (vx_tracer*x_tracer + vy_tracer*y_tracer + vz_tracer*z_tracer)/r_tracer
        
        it_merge = np.argmin( np.abs(t_tracer-t_merge) )
        # it_start = np.argmin( np.abs(vr_tracer-t_merge) )
        # it_diag=len(t_tracer)-1
        it_end = np.argmin( np.abs(t_tracer-0.1) )
        it_diag = np.argmin( np.abs(tem_tracer - tmax_ana) )
        it_dsen_max = 0
        # for it in range(it_diag+1,len(t_tracer)):
        for it in range(it_diag+1,it_end):
            # if sen_tracer[it]>1.0 and dsen_max < np.abs(sen_tracer[it]-sen_tracer[it-1])/sen_tracer[it]:
            if dsen_max < np.abs( (sen_tracer[it]-sen_tracer[it-1])/sen_tracer[it]):
                dsen_max = np.abs( (sen_tracer[it]-sen_tracer[it-1])/sen_tracer[it])
                it_dsen_max = it
                # dsen_max = max([dsen_max,np.abs(sen_tracer[it]-sen_tracer[it-1])/sen_tracer[it]])

        if dsen_max> 0.5 and np.abs((ye_tracer[-1]-ye_tracer[it_dsen_max])/ye_tracer[-1])>0.1:
            mass_ye01_yjump = mass_ye01_yjump + mass_tracer[i]
        
        # it_dsen05=0
        # for it in range(len(t_tracer)-1,0,-1):
        #     if np.abs( (sen_tracer[it]-sen_tracer[it-1])/sen_tracer[it]) > 0.5:
        #         it_dsen05 = it
        #         break
        
        if dsen_max > 0.5:
        #if it_dsen05!=0:
            
            n_tracer = n_tracer + 1

            mass_ye01_sjump = mass_ye01_sjump + mass_tracer[i]

            if np.mod(n_tracer-1,100)==0 :
                color="r"
                alpha=0.5
                zorder=10
                print(ip)
                with open('%06i.dat' % (int(ip)),mode='w') as f:
                    
                    for it in range(len(t_tracer)):
                        
                        ie = max([ied , min([ieu-1, int((np.log10(tem_tracer[it]/11.6e9)-tem_e_min  )*dtemi)])])
                        ie1=ie+1
                        ss = max([0.0, min([1.0 ,       (np.log10(tem_tracer[it]/11.6e9)-tem_e[ie])*dtemi])])
                        ssp= 1.0-ss
                    
                        je = max([jed , min([jeu-1, int((ye_tracer[it]-ye_e_min  )*dyei)])])
                        je1= je +1
                        tt = max([0.0, min([1.0 ,       (ye_tracer[it]-ye_e[je])*dyei])])
                        ttp= 1.0-tt
                        
                        ke = max([ked , min([keu-1, int((np.log10(rho_tracer[it])-rho_e_min  )*drhoi)])])
                        ke1= ke+1
                        uu = max([0.0, min([1.0,        (np.log10(rho_tracer[it])-rho_e[ke])*drhoi])])
                        uup= 1.0-uu
                    
                        zz  = ssp *ttp *uup * zz_e[ke ,je ,ie ]\
                            + ss  *ttp *uup * zz_e[ke ,je ,ie1]\
                            + ssp *tt  *uup * zz_e[ke ,je1,ie ]\
                            + ssp *ttp *uu  * zz_e[ke1,je ,ie ]\
                            + ss  *tt  *uup * zz_e[ke ,je1,ie1]\
                            + ss  *ttp *uu  * zz_e[ke1,je ,ie1]\
                            + ssp *tt  *uu  * zz_e[ke1,je1,ie ]\
                            + ss  *tt  *uu  * zz_e[ke1,je1,ie1]
       
                        aa  = ssp *ttp *uup * aa_e[ke ,je ,ie ]\
                            + ss  *ttp *uup * aa_e[ke ,je ,ie1]\
                            + ssp *tt  *uup * aa_e[ke ,je1,ie ]\
                            + ssp *ttp *uu  * aa_e[ke1,je ,ie ]\
                            + ss  *tt  *uup * aa_e[ke ,je1,ie1]\
                            + ss  *ttp *uu  * aa_e[ke1,je ,ie1]\
                            + ssp *tt  *uu  * aa_e[ke1,je1,ie ]\
                            + ss  *tt  *uu  * aa_e[ke1,je1,ie1]
       
                        xA  = ssp *ttp *uup * xA_e[ke ,je ,ie ]\
                            + ss  *ttp *uup * xA_e[ke ,je ,ie1]\
                            + ssp *tt  *uup * xA_e[ke ,je1,ie ]\
                            + ssp *ttp *uu  * xA_e[ke1,je ,ie ]\
                            + ss  *tt  *uup * xA_e[ke ,je1,ie1]\
                            + ss  *ttp *uu  * xA_e[ke1,je ,ie1]\
                            + ssp *tt  *uu  * xA_e[ke1,je1,ie ]\
                            + ss  *tt  *uu  * xA_e[ke1,je1,ie1]
       
                        f.write( ('%12.4e'*8+'\n') % (t_tracer[it],rho_tracer[it],tem_tracer[it],ye_tracer[it],aa,zz,xA,sen_tracer[it]))

            else:
                color="k"
                alpha=0.1
                zorder=0
                
            ax1.plot(t_tracer-t_merge, r_tracer   ,color=color, ls = "solid", alpha = alpha, linewidth = 1.0, zorder=zorder)
            ax1.scatter(t_tracer[it_dsen_max]-t_merge, r_tracer[it_dsen_max] ,color=color, marker="o", alpha = 0.5, zorder=zorder)
            
            ax2.plot(t_tracer-t_merge, ye_tracer  ,color=color, ls = "solid", alpha = alpha, linewidth = 1.0,zorder=zorder)
            ax2.scatter(t_tracer[it_dsen_max]-t_merge, ye_tracer[it_dsen_max] ,color=color, marker="o", alpha =alpha,zorder=zorder)
            # ax2.plot(t_tracer-t_merge, 4.8e-4/(4e4/vr_tracer) ,color=color, ls = "solid", alpha = alpha, linewidth = 1.0, zorder=zorder)
            
            ax3.plot(t_tracer-t_merge, sen_tracer ,color=color, ls = "solid", alpha = alpha, linewidth = 1.0, zorder=zorder)
            ax3.scatter(t_tracer[it_dsen_max]-t_merge, sen_tracer[it_dsen_max] ,color=color, marker="o", alpha = 0.5, zorder=zorder)
            
            ax4.plot(tem_tracer[it_merge:it_end]/1e9, np.log10(rho_tracer[it_merge:it_end]) ,color=color, ls = "solid", alpha = alpha, linewidth = 1.0, zorder=zorder)
            # ax4.scatter(tem_tracer[it_dsen_max]/1e9, np.log10(rho_tracer[it_dsen_max]) ,color=color, marker="o", alpha = alpha, zorder=zorder)
            ax4.scatter(tem_tracer[it_merge]/1e9, np.log10(rho_tracer[it_merge]) ,color=color, marker="*", alpha = alpha, zorder=zorder)
            
            ax5.plot(tem_tracer[it_merge:it_end]/1e9, sen_tracer[it_merge:it_end] ,color=color, ls = "solid", alpha = alpha, linewidth = 1.0, zorder=zorder)
            # ax5.scatter(tem_tracer[it_dsen_max]/1e9, sen_tracer[it_dsen_max] ,color=color, marker="o", alpha = alpha, zorder=zorder)
            ax5.scatter(tem_tracer[it_merge]/1e9, sen_tracer[it_merge] ,color=color, marker="*", alpha = alpha, zorder=zorder)

            ax6.plot(t_tracer-t_merge, np.log10(rho_tracer) ,color=color, ls = "solid", alpha = alpha, linewidth = 1.0, zorder=zorder)
            ax6.scatter(t_tracer[it_dsen_max]-t_merge, np.log10(rho_tracer[it_dsen_max]) ,color=color, marker="o", alpha = 0.5, zorder=zorder)            
            # ax4.plot(tem_tracer/1e9, ye_tracer ,color="k", ls = "solid", alpha = alpha, linewidth = 1.0)
            # ax4.scatter(tem_tracer[it_dsen_max]/1e9, ye_tracer[it_dsen_max] ,color="k", marker="o", alpha = 0.5)
            # ax4.scatter(tem_tracer[it_merge]/1e9, ye_tracer[it_merge] ,color="k", marker="*", alpha = 0.5)
            
            ax7.plot(t_tracer-t_merge, np.log10(vr_tracer/clight) ,color=color, ls = "solid", alpha = alpha, linewidth = 1.0, zorder=zorder)
            ax7.scatter(t_tracer[it_dsen_max]-t_merge, np.log10(vr_tracer[it_dsen_max]/clight) ,color=color, marker="o", alpha = 0.5, zorder=zorder)
            
    print( ("%6i"+"%12.4e"*2) % (ip,dsen_max,ye_last_tracer[i]))
    

    #ax2.plot(t_tracer-t_merge, hut_tracer ,color="b", ls = "dashed", alpha = 0.1, linewidth = 1.0)
        
print(mass_tracer_total,mass_ye01,mass_ye01_sjump,n_tracer,mass_ye01_yjump)
td=-5e-3;tu=30e-3;scalex="linear"

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
fig1.savefig("tr_sjump_rad.png")

fig2.subplots_adjust(left=0.15,bottom=0.15,right=0.95, top=0.95)
ax2.legend(fontsize=16,frameon=False)
#ax2.set_xlim(td,tu)
ax2.set_xlim(0.001,1)
ax2.set_ylim(0.03,0.1)
#ax2.set_ylim(0.1,300)
ax2.set_xscale('log')
#ax2.set_yscale('log')
ax2.tick_params(axis='both',width=1.5,length=6,which='major',top=True,right=True,pad=10.0)
ax2.tick_params(axis='both',width=1.5,length=3,which='minor',top=True,right=True)
ax2.set_xlabel("$t-t_\\mathrm{merge}$ (s)")
ax2.set_ylabel("$Y_\\mathrm{e}$")
#ax2.set_ylabel("$\\Delta t/(\\Delta x/v^r)$")
#ax2.xaxis.set_major_locator(ticker.MultipleLocator(1.0))
#ax2.xaxis.set_minor_locator(ticker.MultipleLocator(0.2))
fig2.savefig("tr_sjump_ye.png")
#fig2.savefig("tr_sjump_condition.png")

fig3.subplots_adjust(left=0.15,bottom=0.15,right=0.95, top=0.95)
ax3.legend(fontsize=16,frameon=False)
ax3.set_xlim(td,tu)
ax3.set_xscale(scalex)
#ax3.set_ylim(0.01,100)
ax3.set_ylim(-10,30)
#ax3.set_yscale('log')
ax3.tick_params(axis='both',width=1.5,length=6,which='major',top=True,right=True,pad=10.0)
ax3.tick_params(axis='both',width=1.5,length=3,which='minor',top=True,right=True)
ax3.set_xlabel("$t-t_\\mathrm{merge}$ (s)")
ax3.set_ylabel("$s/k_\\mathrm{B}$")
#ax3.xaxis.set_major_locator(ticker.MultipleLocator(1.0))
#ax3.xaxis.set_minor_locator(ticker.MultipleLocator(0.2))
fig3.savefig("tr_sjump_s.png")

fig4.subplots_adjust(left=0.15,bottom=0.15,right=0.95, top=0.95)
ax4.legend(fontsize=16,frameon=False)
ax4.set_xlim(0.3,20)
ax4.set_ylim(3,13)
#ax4.set_xscale(scalex)
ax4.set_xscale('log')
#ax4.set_yscale('log')
ax4.tick_params(axis='both',width=1.5,length=6,which='major',top=True,right=True,pad=10.0)
ax4.tick_params(axis='both',width=1.5,length=3,which='minor',top=True,right=True)
ax4.set_xlabel("$T$ (GK)")
ax4.set_ylabel("Log $\\rho$")
#ax4.xaxis.set_major_locator(ticker.MultipleLocator(2.0))
#ax4.xaxis.set_minor_locator(ticker.MultipleLocator(0.4))
ax4.yaxis.set_major_locator(ticker.MultipleLocator(1.0))
ax4.yaxis.set_minor_locator(ticker.MultipleLocator(0.5))
fig4.savefig("tr_sjump_tem-rho.png")


fig5.subplots_adjust(left=0.15,bottom=0.15,right=0.95, top=0.95)
ax5.legend(fontsize=16,frameon=False)
ax5.set_xlim(0.3,20)
ax5.set_ylim(0.01,100)
#ax5.set_xscale(scalex)
ax5.set_xscale('log')
ax5.set_yscale('log')
ax5.tick_params(axis='both',width=1.5,length=6,which='major',top=True,right=True,pad=10.0)
ax5.tick_params(axis='both',width=1.5,length=3,which='minor',top=True,right=True)
ax5.set_xlabel("$T$ (GK)")
ax5.set_ylabel("$s/k_\\mathrm{B}$")
#ax5.xaxis.set_major_locator(ticker.MultipleLocator(2.0))
#ax5.xaxis.set_minor_locator(ticker.MultipleLocator(0.4))
#ax5.yaxis.set_major_locator(ticker.MultipleLocator(1.0))
#ax5.yaxis.set_minor_locator(ticker.MultipleLocator(0.5))
fig5.savefig("tr_sjump_tem-s.png")


fig6.subplots_adjust(left=0.15,bottom=0.15,right=0.95, top=0.95)
ax6.legend(fontsize=16,frameon=False)
ax6.set_xlim(td,tu)
#ax6.set_ylim(0.01,1)
#ax6.set_xscale(scalex)
#ax6.set_xscale('log')
#ax6.set_yscale('log')
ax6.tick_params(axis='both',width=1.5,length=6,which='major',top=True,right=True,pad=10.0)
ax6.tick_params(axis='both',width=1.5,length=3,which='minor',top=True,right=True)
ax6.set_xlabel("$t$ (s)")
ax6.set_ylabel("$\\rho$")
#ax6.xaxis.set_major_locator(ticker.MultipleLocator(2.0))
#ax6.xaxis.set_minor_locator(ticker.MultipleLocator(0.4))
#ax6.yaxis.set_major_locator(ticker.MultipleLocator(1.0))
#ax6.yaxis.set_minor_locator(ticker.MultipleLocator(0.5))
fig6.savefig("tr_sjump_rho.png")

fig7.subplots_adjust(left=0.15,bottom=0.15,right=0.95, top=0.95)
ax7.legend(fontsize=16,frameon=False)
ax7.set_xlim(td,tu)
ax7.set_ylim(-3,0)
#ax7.set_xscale(scalex)
#ax7.set_xscale('log')
#ax7.set_yscale('log')
ax7.tick_params(axis='both',width=1.5,length=6,which='major',top=True,right=True,pad=10.0)
ax7.tick_params(axis='both',width=1.5,length=3,which='minor',top=True,right=True)
ax7.set_xlabel("$t$ (s)")
ax7.set_ylabel("$v^r/c$")
#ax7.xaxis.set_major_locator(ticker.MultipleLocator(2.0))
#ax7.xaxis.set_minor_locator(ticker.MultipleLocator(0.4))
#ax7.yaxis.set_major_locator(ticker.MultipleLocator(1.0))
#ax7.yaxis.set_minor_locator(ticker.MultipleLocator(0.5))
fig7.savefig("tr_sjump_vel.png")


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

