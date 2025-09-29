import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
#import matplotlib.font_manager as fm
import matplotlib.colors as colors
import matplotlib.cm as cm
import sys
import os
from matplotlib.collections import LineCollection
from mpl_toolkits import mplot3d
import glob
import h5py
from matplotlib.colors import LinearSegmentedColormap

#plt.switch_backend('agg')
#os.environ["PATH"] += os.pathsep + '/data/home/sfujibayashi/libs/texlive/2020/bin/x86_64-linux'


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
    'text.usetex'         : True    }

plt.rcParams.update(params)

cmap_ye = LinearSegmentedColormap.from_list("custom", ["darkred", "red", "orange", "yellow", "green", "white", "blue"], 256)
# cmap_ye= cm.get_cmap('bwr', 16)
norm_ye = colors.Normalize(vmin=0.0, vmax=0.6)

def axisEqual3D(ax):
    extents = np.array([getattr(ax, 'get_{}lim'.format(dim))() for dim in 'xyz'])
    sz = extents[:,1] - extents[:,0]
    centers = np.mean(extents, axis=1)
    maxsize = max(abs(sz))
    r = maxsize/2
    for ctr, dim in zip(centers, 'xyz'):
        getattr(ax, 'set_{}lim'.format(dim))(ctr - r, ctr + r)

len_traj = 50

alpha_grad=np.linspace(0,1,len_traj)

# fn="DD2Tim326_Q4_M135_a75_0056_400m_B5e16_Hon_ng_dyn_unified/traj_00000101.dat"

rlim = 1.0e9

# dir_read = "/scratch/sfujibayashi/Particle_trace_data/data_3D_SFHo_120_150-150mstg_Fugaku/data_sk1_th24_r1.0e+09_test"
dir_read = "/sakura/ptmp/khaya/SFHoTim276_125_165_0025_150mstg_B1e15p05_HLLD_CT_GS-for_PT/Analysis_ptr/data_sk1_th36_r1.0e+09_test"
dir_fig = "./fig"
fn_base = ""

show_time_from_merger = False

if show_time_from_merger:
    t_merge=1.5014e-02
    pass

# fn_list=sorted(glob.glob("DD2Tim326_Q4_M135_a75_0056_400m_B5e16_Hon_ng_dyn_unified/traj_????????.dat"))
fn_list=sorted(glob.glob(dir_read+"/data_???_??????.h5"))

# fn_list = [ "data_sk1_th36_r1.0e+09_test/data_023_000079.h5" ]

# print(len(fn_list))
# sys.exit()

fn = dir_read + "/ana_traj.dat"
ye_particle = np.loadtxt(fn,comments='#',usecols=9)
#print(ye_particle)
#sys.exit()

# fn = "data_post/rejected.dat"
# data=np.loadtxt(fn,comments='#')
# ip_list_rejected=data[:,0]
# #print(len(ip_list_rejected))
# for ip in ip_list_rejected:
#     fn_list.remove("data_post/traj_%08i.dat" % (ip))

#fn_list=[]

# fn_list_dyn=sorted(glob.glob("data_dyn_unified/traj_????????.dat"))

#print(len(fn_list))

# fn = "data_dyn_unified/rejected.dat"
# data=np.loadtxt(fn,comments='#')
# ip_list_rejected=data[:,0]
# #print(len(ip_list_rejected))
# for ip in ip_list_rejected:
#     fn_list_dyn.remove("data_dyn_unified/traj_%08i.dat" % (ip))

# fn_list.extend(fn_list_dyn)

# print("total number of tracers: ",len(fn_list))

fn_list_sub = fn_list[:]
# print(len(fn_list_sub))
n_data = len(fn_list_sub)

skip=1

i_sta = 0
i_offset = 0

# for i, fn in enumerate(fn_list_sub):
for i in range(i_sta,n_data):

    if i%skip!=0:
        continue

    fn = fn_list_sub[i]
    print(i, i+i_offset,fn)
    f0 = h5py.File(fn, 'r')
    time = np.asarray(f0['/time'])[0]

    num_particles = np.asarray(f0['/ipu'])[0]
    x_p  = np.asarray(f0['/x_p'])
    y_p  = np.asarray(f0['/y_p'])
    z_p  = np.asarray(f0['/z_p'])
    print(time,num_particles)
    
    fig = plt.figure(figsize=(10.0, 10.0*0.75))
    # ax = plt.axes(projection="3d", computed_zorder=False)
    ax = plt.axes(projection="3d")
    # ax = fig.add_axes ( (0, 0, 1, 1), projection='3d', computed_zorder=False)

    axisEqual3D(ax)

    x_plot = x_p[::skip]/1e8
    y_plot = y_p[::skip]/1e8
    z_plot = z_p[::skip]/1e8
    r_plot = np.sqrt(x_plot**2 + y_plot**2 + z_plot**2)*1e8
    ye_plot = ye_particle[:num_particles:skip]

    flag_plot = (-rlim/1e8 < x_plot ) & ( x_plot < rlim/1e8 ) & ( -rlim/1e8 < y_plot ) & ( y_plot < rlim/1e8 ) & ( -rlim/1e8 < z_plot ) & ( z_plot <rlim/1e8)
    # flag_plot = (r_plot<rlim*np.sqrt(3.0))
    
    bounds = np.array([rlim/1e8 for ip in range(0,num_particles,skip)])
    
    view = np.array([rlim/1e8, -rlim/1e8, rlim/1e8])
    dist_from_view = np.sqrt( (x_plot-view[0])**2 + (y_plot-view[1])**2 + (z_plot-view[2])**2)
    #print(len(bounds), len(bounds[r_plot<rlim]))
    #sys.exit()
    ax.scatter(x_plot[flag_plot],y_plot[flag_plot],-bounds[flag_plot],color="k", s=1)
    ax.scatter(x_plot[flag_plot], bounds[flag_plot],z_plot[flag_plot],color="k", s=1)
    ax.scatter(-bounds[flag_plot],y_plot[flag_plot],z_plot[flag_plot],color="k", s=1)

    ax.scatter(x_plot[flag_plot], y_plot[flag_plot], z_plot[flag_plot], c=ye_plot[flag_plot], s=40, edgecolor="k", cmap=cmap_ye, norm=norm_ye,alpha=1.0) #,zorder=-dist_from_view)
    
    ax.set_xlim(-rlim/1e8, rlim/1e8)
    ax.set_ylim(-rlim/1e8, rlim/1e8)
    ax.set_zlim(-rlim/1e8, rlim/1e8)
    #ax.set_zlim(-0.0/1e8, rlim/1e8)
    ax.set_xlabel("$x$ (1000 km)", labelpad=15)
    ax.set_ylabel("$y$ (1000 km)", labelpad=15)
    ax.set_zlabel("$z$ (1000 km)", labelpad=15)
    # ax.xaxis.set_major_locator(ticker.MultipleLocator(0.5))
    # ax.xaxis.set_minor_locator(ticker.MultipleLocator(0.1))
    # ax.yaxis.set_major_locator(ticker.MultipleLocator(0.5))
    # ax.yaxis.set_minor_locator(ticker.MultipleLocator(0.1))
    # ax.zaxis.set_major_locator(ticker.MultipleLocator(0.5))
    # ax.zaxis.set_minor_locator(ticker.MultipleLocator(0.1))
    
    #ax.set_zlim(0, np.amax(z_line))
    # sm=cm.ScalarMappable(norm=norm,cmap=cmap)
    # sm.set_array([])
    # cbar=fig.colorbar(sm,ax=ax)
    # cbar.set_label("$Y_\\mathrm{e}$",rotation = 90)

    s_map = cm.ScalarMappable(norm=norm_ye, cmap=cmap_ye)
    s_map.set_array([])
    cbar=fig.colorbar(s_map,pad=0.1, ax=ax)
    #cbar.set_ticks([7.0/6.0, 8.0/6.0, 9.0/6.0])
    #cbar.ax.set_yticklabels(["7/6", "4/3", "3/2"])
    cbar.set_label("$Y_e$",rotation = -90, labelpad= 30)

    if show_time_from_merger:
        ax.text2D(0.05,0.95,'$t-t_\\mathrm{merge}=%7.5f$ s' % (time-t_merge),transform=ax.transAxes,fontsize=20)
    else:
        ax.text2D(0.05,0.95,'$t=%7.5f$ s' % (time),transform=ax.transAxes,fontsize=20)
        pass
    #cbar2.set_clim(0,0.6)
    
    fn = dir_fig + "/%s%06i.png" % (fn_base,i+i_offset)
    fig.savefig(fn)
    plt.close(fig)
    #sys.exit()

sys.exit()
    
