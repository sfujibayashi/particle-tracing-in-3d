import sys
import os
import glob
import h5py

index_loc = 0
index_job_min=1
index_job_max=2
index_jobN_incr=3

# list_locations=[
# ["/scratch/sfujibayashi/DD2Tim326_135_135_0028_12.5mstg_B15.5_HLLD_CT_GS/Analysis_ptr/data_3e8cm", 346, 346, 0], 
# ["/scratch/sfujibayashi/DD2Tim326_135_135_0028_12.5mstg_B15.5_HLLD_CT_GS/Analysis_ptr/data_forward", 347, 347, 0]
# ]
# out_location  = "/scratch/sfujibayashi/DD2Tim326_135_135_0028_12.5mstg_B15.5_HLLD_CT_GS/Analysis_ptr/data_total"

# list_locations=[
# ["/scratch/kiuchi/BHBLpTim326_13625_13625_45km_12.5mstg_B15_HLLD/Analysis_ptr/data_sk1_th36_r1.0e+09_test", 12, 27, 0], 
# ["/scratch/kiuchi/BHBLpTim326_13625_13625_45km_12.5mstg_B15_HLLD_lv16to15/Analysis_ptr/data_sk1_th36_r1.0e+09_test", 28, 56, 0],
# ["/scratch/kiuchi/BHBLpTim326_13625_13625_45km_12.5mstg_B15_HLLD_lv16to15_Mag/Analysis_ptr/data_sk1_th36_r1.0e+09_test", 57, 62, 0],
# ["/scratch/kiuchi/BHBLpTim326_13625_13625_45km_12.5mstg_B15_HLLD_lv15to14_Mag/Analysis_ptr/data_sk1_th36_r1.0e+09_test", 63, 83, 0],
# ["/scratch/kiuchi/BHBLpTim326_13625_13625_45km_12.5mstg_B15_HLLD_lv14to13_Mag_Cowling/Analysis_ptr/data_sk1_th36_r1.0e+09_test", 81, 98, 3]
# ]
# out_location = "/scratch/sfujibayashi/Particle_trace_data/BHBLpTim326_13625_13625_45km_12.5mstg_B15_HLLD/ptr/data_sk1_th36_r1.0e+09_test"


list_locations=[
["/scratch/kiuchi/BHBLpTim326_13625_13625_45km_12.5mstg_B15_HLLD/Analysis_ptr/data_sk1_th36_r3.0e+08_test", 12, 27, 0], 
["/scratch/kiuchi/BHBLpTim326_13625_13625_45km_12.5mstg_B15_HLLD_lv16to15/Analysis_ptr/data_sk1_th36_r3.0e+08_test", 28, 56, 0],
["/scratch/kiuchi/BHBLpTim326_13625_13625_45km_12.5mstg_B15_HLLD_lv16to15_Mag/Analysis_ptr/data_sk1_th36_r3.0e+08_test", 57, 62, 0],
["/scratch/kiuchi/BHBLpTim326_13625_13625_45km_12.5mstg_B15_HLLD_lv15to14_Mag/Analysis_ptr/data_sk1_th36_r3.0e+08_test", 63, 83, 0],
["/scratch/kiuchi/BHBLpTim326_13625_13625_45km_12.5mstg_B15_HLLD_lv14to13_Mag_Cowling/Analysis_ptr/data_sk1_th36_r3.0e+08_test", 81, 98, 3],
["/scratch/sfujibayashi/BHBLpTim326_13625_13625_45km_12.5mstg_B15_HLLD_lv14to13_Mag_Cowling/Analysis_ptr/data_ejecta_forward_3.0e8cm", 99, 103, 3]
]
out_location = "/scratch/sfujibayashi/Particle_trace_data/BHBLpTim326_13625_13625_45km_12.5mstg_B15_HLLD/Analysis_ptr/data_total_3.0e8cm"

n_location=len(list_locations)

# for i in range(n_location):
#     dir_data = list_locations[i][index_loc]
#     cmd = "ls -lh %s/res*h5" % (dir_data)
#     os.system(cmd)
#     pass
# sys.exit()


for i in range(n_location):
    dir_data = list_locations[i][index_loc]
    job_min = list_locations[i][index_job_min]
    job_max = list_locations[i][index_job_max]
    jobN_incr = list_locations[i][index_jobN_incr]
    
    for job in range(job_min, job_max+1,1):
        wildcard_path = "%s/data_%03d_??????.h5" % (dir_data,job)
        list_fn_data = sorted(glob.glob(wildcard_path))

        job_out = job + jobN_incr

        fn_out = out_location + "/steps_%d.dat" % (job_out)
        fn_from = dir_data + "/steps_%d.dat" % (job)
        
        cmd = "ln -s %s %s" % (fn_from, fn_out)
        os.system(cmd)

        ndata = len(list_fn_data)
        print(job, ndata)

        # main body
        for it in range(1, ndata+1):
            print(it)
            fn_out = out_location + "/data_%03d_%06d.h5" % (job_out, it)
            fn_from = dir_data + "/data_%03d_%06d.h5" % (job, it)

            cmd = "ln -s %s %s" % (fn_from, fn_out)
            # print(cmd)
            os.system(cmd)
            # print(fn_out, fn_from)
            # with h5py.File(fn, mode='w') as f:
                
            pass


        #
        
        pass

    pass
