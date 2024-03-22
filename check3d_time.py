import numpy as np
import sys
import h5py
import re

def natural_sort(l): 
    convert = lambda text: int(text) if text.isdigit() else text.lower()
    alphanum_key = lambda key: [convert(c) for c in re.split('([0-9]+)', key)]
    return sorted(l, key=alphanum_key)

#args = sys.argv
#fn = args[1]
clight=2.99792458e10

fn = "/sakura/ptmp/kiuchikn/DD2Tim626_135_135_44km_150mstg_B0_HLLC_FUKA/hdf5_6/vel3d.h5"

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


lv = 14
it = 11647

group = "level%d/data%d" % (lv,it)
print(group)
time = np.asarray(f1[group + "/time"])[0]
vlx = np.asarray(f1[group + "/vx"])/clight
vly = np.asarray(f1[group + "/vy"])/clight
vlz = np.asarray(f1[group + "/vz"])/clight

group = "level%d/data%d" % (lv,it+1)
print(group)
time_prv = np.asarray(f1[group + "/time"])[0]
vlx_b = np.asarray(f1[group + "/vx"])/clight
vly_b = np.asarray(f1[group + "/vy"])/clight
vlz_b = np.asarray(f1[group + "/vz"])/clight

nz = len(vlx)
ny = len(vlx[0])
nx = len(vlx[0][0])
#print(nx,ny,nz)
#sys.exit()

with open("test",mode="w") as f:
    print(time_prv, time)
    l=0
    for k in range(0,ny,2):
        for j in range(0,nx,2):
            f.write( ("%4d"*3+"%12.4E"*6 + "\n") % (j,k,l,vlx[l,k,j],vly[l,k,j],vlz[l,k,j],vlx_b[l,k,j],vly_b[l,k,j],vlz_b[l,k,j]))
sys.exit()
it_lvmin = 1

for i,level in enumerate(list_levels):
    nstep = list_steps[i]
    lv = i+1

    #it = 1
    #!t_start=np.asarray(f1['%s/data%s/time' % (level,it)])[0]
    
    it = nstep
    t_end=np.asarray(f1['%s/data%s/time' % (level,it)])[0]
    
    if lv>lvf:
        it_min = it_lvmin*2**(lv-lvf)
    else:
        it_min = it_lvmin
    #
    
    #print('%s/data%s/time' % (level,it_min))
    t_start=np.asarray(f1['%s/data%s/time' % (level,it_min)])[0]
    
    # take the largest dt and divide by 2 
    #t_end = 0.0
    print("%10s" % (level), "Nstep=%10d,%10d" % (nstep,it_min), "t_start, t_end = %e, %e, (%e)" % (t_start,t_end,t_end-t_start) )
    
    continue
    for it in range(1,nstep+1):
        print(it)
        str1 = '%s/data%d/vx' % (level,it)
        # print(str1)
        vx = np.asarray(f1[str1])
        # print( ("%10d"+"%12.4e"*9) % (it, vx[0,0,0], vx[0,0,-1], vx[0,-1,0], vx[0,-1,-1], vx[-1,-1,0], vx[-1,-1,-1]))
        # print( ("%10d"+"%12.4e") % (it, vx[0,121,121]))

        if vx.any() > 3e10:
            print("some elements > c")
        
# for it in range(1,nstep+1):
#     t=np.asarray(f1['level1/data%s/time' % (it)])[0]
#     print(it,t)
        
