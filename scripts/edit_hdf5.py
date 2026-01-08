import h5py
import numpy as np

fn_hdf = "./res_103.h5"
obj_str = "/scratch/sfujibayashi/BHBLpTim326_13625_13625_45km_12.5mstg_B15_HLLD_lv13_N216_Mag_Cowling/hdf5/103/raw3d.h5"
nstring = len(obj_str)
with h5py.File(fn_hdf, "r+") as f:
    # dset = f["/hdf_read"]
    if "/nstring" in f:
        f["/nstring"][()] = nstring
        pass

    if "/hdf_read" in f:
        del f["/hdf_read"]
        pass

    f.create_dataset(
        "/hdf_read",
        data=np.string_(obj_str),
        dtype=f"S{nstring}"
        #data=obj_str,
        #dtype=h5py.string_dtype(encoding="utf-8")
        #dtype=h5py.special_dtype(vlen=str)
    )
    pass
