import numpy as np
import sys
import h5py
import glob

args = sys.argv
model_dir = args[1]

list_hdf = sorted(glob.glob('%s/hdf5/?'%(model_dir))) + sorted(glob.glob('%s/hdf5/??'%(model_dir))) + sorted(glob.glob('%s/hdf5/???'%(model_dir)))

itt = 0
for dir_hdf in list_hdf:
    fn = dir_hdf + "/raw3d.h5"
    print(fn)
    f1 = h5py.File(fn, 'r')
    list_level=list(f1.keys())
    lv_min=1
    lv_max=len(list_level)
    n_lv=lv_max-lv_min+1
    ft=open("time",mode="w")
    list_obj=list(f1["/level1"].keys())

    nstep=0
    for obj in list_obj:
        if "data" in obj:
            nstep = nstep + 1
    print("Nstep=",nstep)

    t=np.asarray(f1['level1/data%s/time' % (nstep)])[0]
    print(t)
    # for it in range(1,nstep+1):
    #     itt += 1
    #     t=np.asarray(f1['level1/data%s/time' % (it)])[0]
    #     print(it,itt,t)
        
