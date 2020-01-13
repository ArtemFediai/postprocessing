#!/usr/bin/env python3
import argparse
from SimAnalysis import CurrTempSimulation, CurrTempDMRset
import glob
import numpy as np
import matplotlib.pyplot as plt
import os



parser = argparse.ArgumentParser()
parser.add_argument('-marcus', action='store_true', help='Specify if simulation was done with Marcus rates.')
args = parser.parse_args()
print(args)

if args.marcus:
    rates = "Marcus" 
else:
    rates = "Miller"



if not os.path.exists("analysis/Eact_vs_size.txt"):
    sys_dir_unsorted = glob.glob('sys_*')
    list_sys_dirs = sorted(sys_dir_unsorted, key = lambda x: int(x.split('_')[-1]))
    n_sys = len(list_sys_dirs)

    system_sizes = np.zeros(n_sys)
    act_energy   = np.zeros(n_sys)
    Tlow, Thigh  = 250, 350
    for i_s,sys_dir in enumerate(list_sys_dirs):
        simulation = CurrTempSimulation(rates=rates,analyse_Ecoul=False, source_dir=sys_dir+'/',dest_dir='analysis/'+sys_dir)
        
        simulation.get_act_energy(Tlim_low=Tlow,Tlim_high=Thigh)

        simulation.plot_conv_analysis()

        system_sizes[i_s] = simulation.sys_size
        act_energy[i_s]   = simulation.act_energy

    with open("analysis/Eact_vs_size.txt",'w') as f:
        if rates == "Marcus":   
            sim_info = '# Field = {} V/nm, disorder = {} eV, DMR = {} %, lambda = {} eV\n'.format(simulation.field,simulation.dis,simulation.DMR, simulation.lam)
        else:    
            sim_info = '# Field = {} V/nm, disorder = {} eV, DMR = {} %\n'.format(simulation.field,simulation.dis,simulation.DMR)
        f.write(sim_info)
        f.write("# T. range for lin. fit of log J(1000/T): {} - {} K\n".format(Tlow,Thigh))    
        f.write('#\n# size(nm)      Activation energy(meV)\n')
        for i_s in range(n_sys):
            f.write('{0:<8}   {1:<6}\n'.format(system_sizes[i_s],act_energy[i_s]))
    f.close()
else:
    act_energy   = (np.loadtxt("analysis/Eact_vs_size.txt",comments='#',unpack=True))[1]
    system_sizes = (np.loadtxt("analysis/Eact_vs_size.txt",comments='#',unpack=True))[0]
    with open("analysis/Eact_vs_size.txt","r") as f:
        sim_info = (f.readlines())[0]
    sim_info = sim_info.replace('# ','')
plt.plot(system_sizes,act_energy,marker="+", linestyle="None")
plt.xlabel('System size (nm)')
plt.xscale("log")
plt.ylabel('$E_A$ (meV)')
plt.title(sim_info,fontsize=14)
plt.savefig("analysis/Eact_vs_site.png")
plt.close()        
