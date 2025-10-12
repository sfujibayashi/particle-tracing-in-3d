# spript by sho fujibayashi

import os
import subprocess
import glob
import sys

import datetime
import socket

import utils

if_cont=False

# system="yamazaki"; username="sfujibayashi"; work_dir="/scratch/" + username
system="sakura"; username="shofu"; work_dir="/%s/ptmp/" % (system) +username

#model="SFHoTim276_13_14_0025_150mstg_B0_HLLC"; job_min=1; job_max=18; nickname="1314"; coord="STAGGERED"; sym="MIRROR"
#model="SFHoTim276_125_145_0025_200mstg_B0_HLLC"; job_min=12; job_max=12; nickname="125145"; coord="STAGGERED"; sym="MIRROR"
#model="SFHoTim276_125_145_0025_200mstg_B0_HLLC"; job_min=12; job_max=12; nickname="125145l"; coord="STAGGERED"; sym="MIRROR"
#model="SFHoTim276_125_145_0025_250mstg_B0_HLLC"; job_min=14; job_max=14; nickname="125145ll"; coord="STAGGERED"; sym="MIRROR"
# model="DD2Tim326_Q4_M135_a75_0056_270m_B5e16_Hon5png"; job_min=1; job_max=116; nickname="Q4B5H"; coord="NONSTAGGERED"; sym="MIRROR"; dformat="NONFUGAKU"


# # model="DD2Tim326_135_135_0028_12.5mstg_B15.5_HLLD_CT_GS"; job_min=347; job_max=358; nickname="DD2MHD"; coord="STAGGERED"; sym="MIRROR"; dformat="FUGAKU"
# model="DD2Tim326_135_135_0028_12.5mstg_B15.5_HLLD_CT_GS_lv14_to_lv13_Cowling"; job_min=347; job_max=358; nickname="DD2MHD"; coord="STAGGERED"; sym="MIRROR"; dformat="FUGAKU"
# #model="BHBLpTim326_13625_13625_45km_12.5mstg_B15_HLLD"; job_min=88; job_max=88; nickname="BHBLpMHD"; coord="STAGGERED"; sym="MIRROR"; dformat="FUGAKU"

#model="DD2Tim326_135_135_0028_12.5mstg_B15.5_HLLD_CT_GS"; job_min=299; job_max=299; nickname="DD2MHD"; coord="STAGGERED"; sym="MIRROR"; dformat="FUGAKU"
#model="BHBLpTim326_13625_13625_45km_12.5mstg_B15_HLLD_lv14to13_Mag_Cowling"; job_min=88; job_max=88; nickname="BHBLpMHD"; coord="STAGGERED"; sym="MIRROR"; dformat="FUGAKU"
model="DD2Tim326_135_135_0028_12.5mstg_B15.5_HLLD_CT_GS_reduced_lv14_to_lv13_run2"; job_min=344; job_max=346; nickname="DD2MHD"; coord="STAGGERED"; sym="MIRROR"; dformat="FUGAKU"

fn_eos = "/sakura/ptmp/shofu/EOS/EOS_Hempel_DD2Tim326_TF"; nrho=426; nye=60; ntemp=131
#fn_eos = "/sakura/ptmp/shofu/EOS/EOS_Hempel_SFHoTim326_TF"; nrho=408; nye=60; ntemp=131
#fn_eos = "/scratch/sfujibayashi/EOS/EOS_Hempel_DD2Tim_TF_326"; nrho=426; nye=60; ntemp=131
#fn_eos = "/scratch/sfujibayashi/EOS/EOS_BHBLpTim_rho453_temp156_ye061_ierd076_knuc376"; nrho=453; nye=60; ntemp=156

#dir_read = work_dir + "/" + model + "/hdf5"
#dir_read = "/sakura/ptmp/khaya/" + model + "/hdf5/"
dir_read = "/sakura/ptmp/kiuchikn/" + model + "/hdf5_"
#dir_read = "/scratch/kiuchi/BHBLpTim326_13625_13625_45km_12.5mstg_B15_HLLD_lv14to13_Mag_Cowling/hdf5/"

# parameters #
n_theta = 9
#n_theta = 64
it_start= 0
it_skip = 1
it_skip_out = 1 #it_skip

rfl=3e8
rin=0.0
mass_crit = 1e-5
mass_min  = 1e-12
backward="T"
volumebased="F"

restart="F"

incr_next=10

info="3e8km"
#info="sk%i_th%i_r%7.1e_omp" % (it_skip,n_theta,rfl)
#info="bind_%ims%i" % (it_skip,n_theta)
#info = "forward_close"
#info = "ana"

print("sub info = ",info)
if info!="":
   info = "_"+info

with open("src/macro.h",mode="w") as f:
   f.write("#define %s\n" % (coord))
   f.write("#define %s\n" % (sym))
   f.write("#define %s\n" % (dformat))
   pass

if not os.path.exists(fn_eos):
    print("EOS file does not exists!")
    sys.exit()

queue="p.sakura"
misc=""
nodes=1
MPI=1
MPI_per_node=1
CPUs_per_task=40
OMP_NUM_THREADS=CPUs_per_task
whour=24
wmin=00
wsecond=00
compf="h5pfc"
option="-convert big_endian -mcmodel=large -shared-intel -fpic -qopenmp -xCORE-AVX512 -qopt-zmm-usage=high"

prog="ptr.out"

cwdir = os.getcwd()
current_time = datetime.datetime.now()
# os_name = os.name
host_name = socket.gethostname()

stamp = "# Generated for model: %s" % (model) + ", on: %s" % (current_time) + ", in the host: %s." % (host_name)


with  open("make.inc", mode="w") as f_make_inc:
   f_make_inc.write("#%s\n" % (stamp))
   f_make_inc.write("COMPILER_NAME = %s\n" % (compf))
   f_make_inc.write("OPTION = %s\n" % (option))
   pass

# fn_makefile = "makefile_raw"
# with open(fn_makefile,mode="r") as f:
#     lines = f.readlines()

# with open("makefile",mode="w") as f:
#     for line in lines:
#         line_rep = line.replace("XXXXX",prog).replace("YYYYY",option).replace("COMPF",compf).replace("COMPC",compc)
#         f.write(line_rep)

print( "queue :", queue)
print( "model :", model)
print( "job    : %i - %i" % (job_min, job_max) )
print( "wall time : %i:%i:%i" % (whour,wmin,wsecond) )
print( "nodes = %i" % (nodes) )
print( "MPI processes per node :",MPI_per_node)
print( "CPU per task           :",CPUs_per_task)

dir_job = work_dir + "/" + model + "/Analysis_ptr"
if not os.path.exists(dir_job):
    cmd = "mkdir -p %s" % (dir_job)
    os.system(cmd)
    pass

dir_err = work_dir + "/" + model + "/Analysis_ptr/err"
if not os.path.exists(dir_err):
    cmd = "mkdir %s" % (dir_err)
    os.system(cmd)
    pass

dir_std = work_dir + "/" + model + "/Analysis_ptr/out"
if not os.path.exists(dir_std):
   cmd = "mkdir %s" % (dir_std)
   os.system(cmd)
   pass

dir_exists=False
dir_out  = work_dir + "/" + model + "/Analysis_ptr/data" + info
# dir_out  = "./"
if os.path.exists(dir_out):
    print("directory exists :",dir_out)
    dir_exists=True
else:
    cmd = "mkdir %s" % (dir_out)
    os.system(cmd)
    pass

if dir_exists:
    run=input('Directory already exists. Overwrite parameter files? [y/n] :')
    if run=='n' or run=='N':
        sys.exit()
    elif run=='y' or run=='Y':
       print("will be overwritten.")
    else:
       print("not specified. stop.")
       pass
    pass

# with open(dir_job + "/para.dat",mode="w") as f:
#     f.write("# result saved in:\n")
#     f.write("%s\n" % (dir_out))

with open(dir_out + "/parameters.dat",mode="w") as f:
    f.write("# model name:\n")
    f.write("%s\n" % (model))
    f.write("# job_min:\n")
    f.write("%i\n" % (job_min))
    f.write("# job_max:\n")
    f.write("%i\n" % (job_max))
    
    f.write("# data read from:\n")
    f.write("%s\n" % (dir_read))
    f.write("# result saved in:\n")
    f.write("%s\n" % ("./"))
    
    f.write("# Evolving backward?:\n")
    f.write("%s\n" % (backward))
    f.write("# Particles are set in a volume-based way?:\n")
    f.write("%s\n" % (volumebased))

    f.write("# number of snapshort where the trace starts (from the first/last if it is 0):\n")
    f.write("%i\n" % (it_start))
    f.write("# interval of trace:\n")
    f.write("%i\n" % (it_skip))
    f.write("# interval of output:\n")
    f.write("%i\n" % (it_skip_out))
    f.write("# angular resolution for particle setting:\n")
    f.write("%i\n" % (n_theta))
    f.write("# R_out (where flux-based particles are set):\n")
    f.write("%e\n" % (rfl))
    f.write("# R_in:\n")
    f.write("%e\n" % (rin))
    f.write("# M_crit (volume-based particle mass should be below):\n")
    f.write("%e\n" % (mass_crit))
    f.write("# M_min  (volume-based particle mass should be above):\n")
    f.write("%e\n" % (mass_min))

    f.write("# EOS file:\n")
    f.write("%s\n" % (fn_eos))
    f.write("# index of rho, ye, temp:\n")
    f.write("%i %i %i\n" % (nrho,nye,ntemp))
    f.write("# incrementation in the next job:\n")
    f.write("%i\n" % (incr_next))

    f.write("# Restart flag:\n")
    f.write("%s\n" % (restart))
    f.write("# job number of checkpoint file:\n")
    f.write("%i\n" % (0))
    f.write("# time step of checkpoint file:\n")
    f.write("%i\n" % (0))
    pass

with open(dir_out + "/ptr.para",mode="w") as f:
    f.write("# model name:\n")
    f.write("model = \"%s\"\n" % (model))
    f.write("# job_min:\n")
    f.write("job_min = %d\n" % (job_min))
    f.write("# job_max:\n")
    f.write("job_max = %d\n" % (job_max))
    
    f.write("# data read from:\n")
    f.write("dir_read = \"%s\"\n" % (dir_read))
    f.write("# result saved in:\n")
    f.write("dir_out  = \"%s\"\n" % (dir_out))
    
    f.write("# Evolving backward?:\n")
    f.write("mode_backward = %s\n" % (backward))
    f.write("# Particles are set in a volume-based way?:\n")
    f.write("mode_volbased = %s\n" % (volumebased))

    f.write("# number of snapshort where the trace starts (from the first/last if it is 0):\n")
    f.write("it_start = %d\n" % (it_start))
    f.write("# interval of trace:\n")
    f.write("it_skip = %d\n" % (it_skip))
    f.write("# interval of output:\n")
    f.write("it_skip_out = %d\n" % (it_skip_out))
    f.write("# angular resolution for particle setting:\n")
    f.write("n_theta = %d\n" % (n_theta))
    f.write("# R_out (where flux-based particles are set):\n")
    f.write("rfl = %e\n" % (rfl))
    f.write("# R_in:\n")
    f.write("rin = %e\n" % (rin))
    f.write("# M_crit (volume-based particle mass should be below):\n")
    f.write("mass_crit = %e\n" % (mass_crit))
    f.write("# M_min  (volume-based particle mass should be above):\n")
    f.write("mass_min = %e\n" % (mass_min))

    f.write("# EOS file:\n")
    f.write("fn_eos = \"%s\"\n" % (fn_eos))
    f.write("# index of rho, ye, temp:\n")
    f.write("nrho = %d\n" % (nrho))
    f.write("ntemp= %d\n" % (ntemp))
    f.write("nye  = %d\n" % (nye))
    f.write("# incrementation in the next job:\n")
    f.write("incr_next = %d\n" % (incr_next))

    f.write("# Restart flag:\n")
    f.write("restart = %s\n" % (restart))
    f.write("# job number of checkpoint file:\n")
    f.write("job_restart = %d\n" % (0))
    f.write("# time step of checkpoint file:\n")
    f.write("it_restart = %d\n" % (0))
    pass

    
with open("sub_script_sakura",mode="r") as f:
    lines=f.readlines()
    pass

list_replace=[
["${dir_std}",dir_std],
["${dir_err}",dir_err],
["${nickname}",nickname],
["${cjob}","ptr"],
["${queue}",queue],
["${nodes}","%i" % (nodes)],
["${MPI_per_node}","%i" % (MPI_per_node)],
["${CPUs_per_task}","%i" % (CPUs_per_task)],
["${mail}","sho.fujibayashi@gmail.com"],
["${hour}","%02i" % (whour)],
["${min}","%02i" % (wmin)],
["${second}","%02i" % (wsecond)],
["${second}","%02i" % (wsecond)],
["${dir_out}",dir_out],
["${work_dir}",work_dir],
["${prog}",prog],
["${exe}","ptr.out"],
["${OMP_NUM_THREADS}","%i" % (OMP_NUM_THREADS)],
["${model}",model],
["${info}",info],
["${misc}",misc]
]


#list_job = [job for job in range(job_min,job_max+1)]
njob = int((job_max-job_min+1)/incr_next) + 1


for job in range(1,njob+1):
   with open("sub_ptr_%s%s_%i.sh" % (model,info,job) ,mode="w") as f:
      for line in lines:
         line1=line
         for rep_pair in list_replace:
            #print(rep_pair[0],rep_pair[1])
            line1=line1.replace(rep_pair[0],rep_pair[1])
         line1=line1.replace("${job}","%i" %(job))
         f.write(line1)
         pass
      pass
   pass
# Ask compile or not

if_compile=input('Compile the code? [y/n] :')
if if_compile=='n' or if_compile=='N':
   print("will not be compiled.")
elif if_compile=='y' or if_compile=='Y':
   
   cmd="make clean"
   stdout, stderr, return_code = utils.exec_subprocess(cmd)
   
   print("Compiling...")
   cmd="make"
   stdout, stderr, return_code = utils.exec_subprocess(cmd)
   print(stdout.decode())
   print(stderr.decode())

   if return_code!=0:
      print("compile failed. stop.")
      sys.exit()
      pass
else:
   print("not specified. stop.")
   pass

# Ask move the exec or not

if if_compile=='y' or if_compile=='Y':
   run=input('move the executables to the output directory? [y/n] :')
else:
   run="n"
if run=='n' or run=='N':
   print("exec. will not be moved.")
elif run=='y' or run=='Y':

   cmd="mv ./bin/%s %s/" %( prog, dir_out )
   stdout, stderr, return_code = utils.exec_subprocess(cmd)
   print(stdout.decode())
   print(stderr.decode())

   if return_code!=0:
      print("command failed. stop.")
      sys.exit()
      pass
else:
   print("not specified. stop.")
   pass
   
# Ask submit or not

run=input('Submit job? [y/n] :')
if run=='n' or run=='N':
   print("will not be submitted")
elif run=='y' or run=='Y':
   print(njob)
   #for job in list_job:
   for job in range(1,njob+1):
      
      fn_sub = "sub_ptr_%s%s_%i.sh" % (model,info,job)
      print(job, fn_sub)
      if if_cont or job>1:
         pid = utils.find_id(username, nickname, job-1)
         if pid==0:
            print("running job not found. stop")
            sys.exit()
         else:
            print("job = ",job)
            cmd = "sbatch -d, --dependency=afterok:%i %s" % (pid,fn_sub)
            os.system(cmd)
      else:
         cmd = "sbatch %s" % (fn_sub)
         print(cmd)
         os.system(cmd)
         pass
      pass
else:
   print("not specified. stop.")
   sys.exit()
   
