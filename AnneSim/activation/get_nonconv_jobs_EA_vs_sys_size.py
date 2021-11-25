#!/usr/bin/env python3
from SimAnalysis import CurrTempSimulation, CurrTempDMRset
import glob
import numpy as np
import argparse
import os
import shutil

parser = argparse.ArgumentParser()
parser.add_argument('-marcus', action='store_true', help='Specify if simulation was done with Marcus rates.')
args = parser.parse_args()
print(args)

if args.marcus:
    rates = "Marcus" 
else:
    rates = "Miller"

sys_dir_unsorted = glob.glob('sys_*')
sys_dirs = sorted(sys_dir_unsorted, key = lambda x: int(x.split('_')[-1]))

for i_s, sys_dir in enumerate(sys_dirs):
    sim = CurrTempSimulation(rates=rates,analyse_Ecoul=False, source_dir=sys_dir+'/',dest_dir='analysis/'+sys_dir)
        
    converged_jobs = np.loadtxt("analysis/"+sys_dir+"/conv_analysis.txt",unpack=True)[3]
    nonconv_jobs = []
    n_temp  = len(sim.temp)
    for i_t in range(n_temp):
        if converged_jobs[i_t] != sim.n_r:
            n_cpus = int(np.loadtxt("joblist_sys_"+str(i_s),usecols=0)[0])
            run_kmc_name = np.loadtxt("joblist_sys_"+str(i_s),usecols=1,dtype=str)[0]
            for i_r in range(sim.n_r):
                nonconv_jobs.append("{} {} {} {} {}\n".format(n_cpus, run_kmc_name, i_s, i_t, i_r))
            print("sys_{}/temp_{} not fully converged. Moving files to non_conv".format(i_s,i_t)) 
            if not os.path.exists("non_conv"):
                os.makedirs("non_conv")   
            try:
                shutil.move("sys_{}/temp_{}".format(i_s,i_t),"non_conv/sys_{}/temp_{}".format(i_s,i_t))
            except:
                print("Already moved")
                    
        if len(nonconv_jobs) != 0:
            print("new_joblist_sys_{} and new .sh file are being generated.".format(i_s))

            new_joblist = open("new_joblist_sys_"+str(i_s),'w')
            for line in nonconv_jobs: new_joblist.write(line)
            new_joblist.close()

            old_sh_name = glob.glob("submit_*_sys_{}.sh".format(i_s))[0]
            old_sh_file = open(old_sh_name,'r')
            sh_data = old_sh_file.readlines()
            old_sh_file.close()

            sh_data[-1] = sh_data[-1].replace("joblist","new_joblist")
            new_sh_name = old_sh_name.replace("submit","new_submit")
            new_sh_file = open(new_sh_name,"w")
            for line in sh_data:
                new_sh_file.write(line)
            new_sh_file.close()