import numpy as np
import matplotlib.pyplot as plt
#import matplotlib.ticker as ticker
#import matplotlib.font_manager as fm
import matplotlib.colors as colors
import matplotlib.cm as cm
import sys
#import os
from matplotlib.collections import LineCollection
from mpl_toolkits import mplot3d
import glob


len_traj = 50

alpha_grad=np.linspace(0,1,len_traj)

# fn="DD2Tim326_Q4_M135_a75_0056_400m_B5e16_Hon_ng_dyn_unified/traj_00000101.dat"

fn_list=sorted(glob.glob("DD2Tim326_Q4_M135_a75_0056_400m_B5e16_Hon_ng_dyn_unified/traj_????????.dat"))
# fn_list=[
# "DD2Tim326_Q4_M135_a75_0056_400m_B5e16_Hon_ng_dyn_unified/traj_00000001.dat",
# "DD2Tim326_Q4_M135_a75_0056_400m_B5e16_Hon_ng_dyn_unified/traj_00000101.dat",
# "DD2Tim326_Q4_M135_a75_0056_400m_B5e16_Hon_ng_dyn_unified/traj_00000201.dat",
# "DD2Tim326_Q4_M135_a75_0056_400m_B5e16_Hon_ng_dyn_unified/traj_00000301.dat",
# "DD2Tim326_Q4_M135_a75_0056_400m_B5e16_Hon_ng_dyn_unified/traj_00000401.dat",
# "DD2Tim326_Q4_M135_a75_0056_400m_B5e16_Hon_ng_dyn_unified/traj_00000501.dat",
# "DD2Tim326_Q4_M135_a75_0056_400m_B5e16_Hon_ng_dyn_unified/traj_00000601.dat"]


t_line=[]
x_line=[]
y_line=[]
z_line=[]
ye_line=[]
for fn in fn_list:
    print(fn)
    data=np.loadtxt(fn,comments='#')
    

    t_line.append ( data[:,0])
    x_line.append ( data[:,1])
    y_line.append ( data[:,2])
    z_line.append ( data[:,3])
    ye_line.append( data[:,9])
    print(len(t_line[-1]))
    # ye_fin=ye_line[-1]



norm = colors.Normalize(vmin=0,vmax=0.6)
cmap = plt.get_cmap('jet_r')

#print(cm.jet_r(ye_fin/0.6)[:3])
#sys.exit()


i_plot = 0
for it in range(len_traj,len(t_line[0]),3):

    fig = plt.figure()
    ax = plt.axes(projection="3d")
    
    for ip in range(len(t_line)):

        if it< len(t_line[ip]):
            points3D=np.array([x_line[ip][it-len_traj:it],y_line[ip][it-len_traj:it],z_line[ip][it-len_traj:it]]).T.reshape(-1,1,3)
            segments=np.concatenate([points3D[:-1], points3D[1:]], axis=1)
            
            alphas=np.linspace(0,1,len(segments))
            
            rgba_colors = np.zeros((len(segments), 4))
            rgba_colors[:,0:3] = cm.jet_r(ye_line[ip][it]/0.6)[0:3]#[1.0,0.0,0.0]
            rgba_colors[:,3] = alphas
            
            lc = mplot3d.art3d.Line3DCollection(segments, colors=rgba_colors)
            #sys.exit()
            #fig,a = plt.subplots()
            ax.add_collection(lc)


    # # len_traj=10
    # it_plot_d = max([0,it-len_traj])
    # it_plot_u = it
    # itt=it
    # ax.plot3D(x_line[itt:itt+len_traj], y_line[itt:itt+len_traj], z_line[itt:itt+len_traj], 'k',alpha=alpha_grad)
    # # for itt in range(it_plot_d,it_plot_u):
    # #     alpha=float(itt-it_plot_d)/float(it_plot_u-it_plot_d)
    # #     #print(itt)
    # #     ax.plot3D(x_line[itt:itt+2], y_line[itt:itt+2], z_line[itt:itt+2], 'k',alpha=alpha)
    # #ax.plot3D(x_line[itt:itt+len_traj]-0.1, y_line[itt:itt+len_traj], z_line[itt:itt+len_traj], 'k',alpha=0.3)

    ax.set_xlim(-1e9, 1e9)
    ax.set_ylim(-1e9, 1e9)
    ax.set_zlim(0, 1e9)
    #ax.set_zlim(0, np.amax(z_line))
    cbar2=fig.colorbar(cm.ScalarMappable(norm=norm,cmap=cmap),ax=ax)
    #cbar2.set_clim(0,0.6)
    
    fig.savefig("%05i.png"%(i_plot))
    plt.close(fig)
    print(it,len(t_line))
    i_plot = i_plot + 1

    #break

plt.show()
