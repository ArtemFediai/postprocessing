#!/usr/bin/env python3
import argparse
from SimAnalysis import CurrTempSimulation, CurrTempDMRset
import glob
import numpy as np
import matplotlib.pyplot as plt
import os
import shutil



parser = argparse.ArgumentParser()
parser.add_argument('-ect', action='store_true', help='Specify if EA vs. Ect should be determined.')
parser.add_argument('-dis', action='store_true', help='Specify if EA vs. disorder should be determined.')
parser.add_argument('-lam', action='store_true', help='Specify if EA vs. lambda should be determined.')

parser.add_argument('-nonconv', action='store_true', help='get non conv. jobs.')
args = parser.parse_args()
print(args)

rates = "Marcus" 

if args.ect:
    dep = "Ecoul"
elif args.dis:
    dep = "dis"
elif args.lam:
    dep = "l"

if not args.nonconv:
    if not os.path.exists("analysis/Eact_vs_{}.txt".format(dep)):
        dir_unsorted = glob.glob('{}_*'.format(dep))
        list_dirs = sorted(dir_unsorted, key = lambda x: int(x.split('_')[-1]))
        n_x = len(list_dirs)

        x_value      = np.zeros(n_x)
        act_energy   = np.zeros(n_x)
        Tlow, Thigh  = 250, 350
        for i, xdir in enumerate(list_dirs):
            simulation = CurrTempSimulation(rates=rates,analyse_Ecoul=args.ect, source_dir=xdir+'/',dest_dir='analysis/'+xdir)
            
            simulation.get_act_energy(Tlim_low=Tlow,Tlim_high=Thigh)

            simulation.plot_conv_analysis()

            if args.ect:
                x_value[i]      = simulation.Ecoul
            elif args.dis:
                x_value[i]      = simulation.dis
            elif args.lam:
                x_value[i]      = simulation.lam        
            
            act_energy[i]   = simulation.act_energy

        with open("analysis/Eact_vs_{}.txt".format(dep),'w') as f:    
            f.write("# T. range for lin. fit of log J(1000/T): {} - {} K\n".format(Tlow,Thigh))
            if args.ect:
                f.write('#\n# Ect(eV)      Activation energy(meV)\n')
            elif args.dis:
                f.write('#\n# Disorder(eV)      Activation energy(meV)\n')
            elif args.lam:
                f.write('#\n# lambda(eV)      Activation energy(meV)\n')
            for i in range(n_x):
                f.write('{0:<8}   {1:<6}\n'.format(x_value[i],act_energy[i]))
        f.close()
    else:
        act_energy = (np.loadtxt("analysis/Eact_vs_{}.txt".format(dep),comments='#',unpack=True))[1]
        x_value    = (np.loadtxt("analysis/Eact_vs_{}.txt".format(dep),comments='#',unpack=True))[0]
    plt.plot(x_value,act_energy,marker="+", linestyle="None")
    if args.ect:
        plt.xlabel('$E_C$ (eV)')
    elif args.dis:
        plt.xlabel('$\sigma$ (eV)')
    elif args.lam:
        plt.xlabel('$\lambda$ (eV)')
    plt.ylabel('$E_A$ (meV)')
    plt.savefig("analysis/Eact_vs_{}.png".format(dep))
    plt.close()        
else:
    dir_unsorted = glob.glob('{}_*'.format(dep))
    list_dirs = sorted(dir_unsorted, key = lambda x: int(x.split('_')[-1]))     
    for i, xdir in enumerate(list_dirs):
        converged_jobs = np.loadtxt("analysis/"+xdir+"/conv_analysis.txt",unpack=True)[3]
        nonconv_jobs = []
        simulation = CurrTempSimulation(rates=rates,analyse_Ecoul=args.ect, source_dir=xdir+'/',dest_dir='analysis/'+xdir)
        n_t = len(simulation.temp)
        for i_t in range(n_t):
            if converged_jobs[i_t] != simulation.n_r:
                n_cpus = int(np.loadtxt("joblist_{}_{}".format(dep,i),usecols=0)[0])
                run_kmc_name = np.loadtxt("joblist_{}_{}".format(dep,i),usecols=1,dtype=str)[0]
                for i_r in range(simulation.n_r):
                    nonconv_jobs.append("{} {} {} {} {}\n".format(n_cpus, run_kmc_name, i, i_t, i_r))
                print("{}_{}/temp_{} not fully converged. Moving files to non_conv".format(dep,i,i_t)) 
                if not os.path.exists("non_conv"):
                    os.makedirs("non_conv")   
                try:
                    shutil.move("{}_{}/temp_{}".format(dep,i,i_t),"non_conv/{}_{}/temp_{}".format(dep,i,i_t))
                except:
                    print("Already moved")
                
        if len(nonconv_jobs) != 0:
            print("new_joblist_{}_{} and new .sh file are being generated.".format(dep,i))

            new_joblist = open("new_joblist_{}_{}".format(dep,i),'w')
            for line in nonconv_jobs: new_joblist.write(line)
            new_joblist.close()

            old_sh_name = glob.glob("submit_EA_*{}_{}.sh".format(dep,i))[0]
            old_sh_file = open(old_sh_name,'r')
            sh_data = old_sh_file.readlines()
            old_sh_file.close()

            sh_data[-1] = sh_data[-1].replace("joblist","new_joblist")
            new_sh_name = old_sh_name.replace("submit","new_submit")
            new_sh_file = open(new_sh_name,"w")
            for line in sh_data:
                new_sh_file.write(line)
            new_sh_file.close()