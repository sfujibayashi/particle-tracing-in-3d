import os
import subprocess
import glob
import sys

def exec_subprocess(cmd: str, raise_error=True):
   
   child = subprocess.Popen(cmd, shell=True,
                            stdout=subprocess.PIPE, stderr=subprocess.PIPE)
   stdout, stderr = child.communicate()
   return_code = child.returncode
   if return_code != 0 and raise_error:
      raise Exception(f"command return code is not 0. got {return_code}. stderr = {stderr}")
   
   return stdout, stderr, return_code


def find_id(username: str, nickname: str, job: int):
   cmd = "squeue -u %s -o \"%%7A %%20u %%30j\"" % (username)
   #print(cmd)
   child = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
   stdout,stderr = child.communicate()
   return_code = child.returncode
   
   pid = 0
   
   str1=stdout.decode()
   lines=str1.split("\n")
   pname = nickname+"-%i" % (job)
   for line in lines:
      if pname in line:
         pid = int(line.split()[0])
         break
   return pid
