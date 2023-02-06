import os
import subprocess
import glob
import sys

# spript by sho
import utils

if_cont=False

model="SFHoTim276_13_14_0025_150mstg_B0_HLLC"; job_min=6; job_max=6; nickname="1314"; coord="STAGGERED"; sym="MIRROR"
#model="SFHoTim276_125_145_0025_200mstg_B0_HLLC"; job_min=12; job_max=12; nickname="125145"; coord="STAGGERED"; sym="MIRROR"
#model="SFHoTim276_125_145_0025_200mstg_B0_HLLC"; job_min=12; job_max=12; nickname="125145l"; coord="STAGGERED"; sym="MIRROR"
#model="SFHoTim276_125_145_0025_250mstg_B0_HLLC"; job_min=14; job_max=14; nickname="125145ll"; coord="STAGGERED"; sym="MIRROR"

# parameters #
n_theta = 32
#n_theta = 64
it_start= 18
it_skip = 1
it_skip_out = 1 #it_skip
rfl=5e8
rin=1e7
mass_crit = 1e-6
mass_min  = 1e-12
backward="T"

restart="N"

#info="sk%i_th%i_r%7.1e_omp" % (it_skip,n_theta,rfl)
#info="bind_%ims%i" % (it_skip,n_theta)
#info = "forward_close"
info = "ana"

print(info)

#fn_eos = "/sakura/ptmp/shofu/EOS/EOS_Hempel_DD2Tim326_TF"; nrho=426; nye=60; ntemp=131
fn_eos = "/sakura/ptmp/shofu/EOS/EOS_Hempel_SFHoTim326_TF"; nrho=408; nye=60; ntemp=131


with open("macro.h",mode="w") as f:
   f.write("#define %s\n" % (coord))
   f.write("#define %s\n" % (sym))

if not os.path.exists(fn_eos):
    print("EOS file does not exists!")
    sys.exit()

machine="sakura"
username="shofu"
work_dir="/sakura/ptmp/"+username

queue="p.sakura"
nodes=1
MPI=1
MPI_per_node=1
CPUs_per_task=40
OMP_NUM_THREADS=CPUs_per_task
whour=24
wmin=00
wsecond=00
compf="h5pfc"
compc="h5pcc"
option="-convert big_endian -mcmodel=large -shared-intel -fpic -qopenmp -xCORE-AVX512 -qopt-zmm-usage=high"

prog="ptr_"+model+".out"

fn_makefile = "makefile_raw"
with open(fn_makefile,mode="r") as f:
    lines = f.readlines()

with open("makefile",mode="w") as f:
    for line in lines:
        line_rep = line.replace("XXXXX",prog).replace("YYYYY",option).replace("COMPF",compf).replace("COMPC",compc)
        f.write(line_rep)


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

dir_err = work_dir + "/" + model + "/Analysis_ptr/err"
if not os.path.exists(dir_err):
    cmd = "mkdir %s" % (dir_err)
    os.system(cmd)
dir_std = work_dir + "/" + model + "/Analysis_ptr/out"
if not os.path.exists(dir_std):
    cmd = "mkdir %s" % (dir_std)
    os.system(cmd)

dir_exists=False
dir_out  = work_dir + "/" + model + "/Analysis_ptr/data_" + info
if os.path.exists(dir_out):
    print("directory exists :",dir_out)
    dir_exists=True
else:
    cmd = "mkdir %s" % (dir_out)
    os.system(cmd)

dir_read = work_dir + "/" + model + "/hdf5"
#dir_read = "/sakura/ptmp/kiuchikn/" + model + "/hdf5"

if dir_exists:
    run=input('Directory already exists. Overwrite parameter files? [y/n] :')
    if run=='n' or run=='N':
        sys.exit()
    elif run=='y' or run=='Y':
       print("will be overwritten.")
    else:
       print("not specified. stop.")

with open(dir_job + "/para.dat",mode="w") as f:
    f.write("# result saved in:\n")
    f.write("%s\n" % (dir_out))

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
    f.write("%s\n" % (dir_out))
    
    f.write("# Evolving backward?:\n")
    f.write("%s\n" % (backward))

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

    f.write("# Restart flag:\n")
    f.write("%s\n" % (restart))
    f.write("# job number of checkpoint file:\n")
    f.write("%i\n" % (0))
    f.write("# time step of checkpoint file:\n")
    f.write("%i\n" % (0))


# with open(dir_out + "/restart_info.dat",mode="w") as f:
#     f.write("%i\n" % (job_max))
#     f.write("%i\n" % (job_max))
#     f.write("N\n")
    
with open("sub_script_sakura",mode="r") as f:
    lines=f.readlines()
    

list_replace=[
["${dir_out}",dir_std],
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
["${dir_job}",dir_job],
["${work_dir}",work_dir],
["${prog}",prog],
["${exe}","a2.out"],
["${OMP_NUM_THREADS}","%i" % (OMP_NUM_THREADS)],
["${model}",model],
["${info}",info]
]


list_job = [job for job in range(job_min,job_max+1)]
njob = job_max-job_min+1

for job in list_job:
   with open("sub_ptr_%s_%s_%i.sh" % (model,info,job) ,mode="w") as f:
      for line in lines:
         line1=line
         for rep_pair in list_replace:
            #print(rep_pair[0],rep_pair[1])
            line1=line1.replace(rep_pair[0],rep_pair[1])
         line1=line1.replace("${job}","%i" %(job))
         f.write(line1)
         

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

else:
   print("not specified. stop.")


# Ask move the exec or not

if if_compile=='y' or if_compile=='Y':
   run=input('move the executables to bin? [y/n] :')
else:
   run="n"
if run=='n' or run=='N':
   print("exec. will not be moved.")
elif run=='y' or run=='Y':

   cmd="mv %s %s/bin/" %( prog, work_dir )
   stdout, stderr, return_code = utils.exec_subprocess(cmd)
   print(stdout.decode())
   print(stderr.decode())

   if return_code!=0:
      print("command failed. stop.")
      sys.exit()

else:
   print("not specified. stop.")

   
# Ask submit or not

run=input('Submit job? [y/n] :')
if run=='n' or run=='N':
   print("will not be submitted")
elif run=='y' or run=='Y':
   #for job in list_job:
   for job in range(job_max,job_min-1,-1):
      
      fn_sub = "sub_ptr_%s_%s_%i.sh" % (model,info,job)

      if if_cont or job!=job_max:
         pid = utils.find_id(username, nickname, job+1)
         if pid==0:
            print("running job not found. stop")
            sys.exit()
         else:
            print("job = ",job)
            cmd = "sbatch -d, --dependency=afterok:%i %s" % (pid,fn_sub)
            os.system(cmd)
      else:
         cmd = "sbatch %s" % (fn_sub)
         os.system(cmd)
         
else:
   print("not specified. stop.")
   sys.exit()
