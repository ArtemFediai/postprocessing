#!/usr/bin/env python3
import numpy as np 
import glob
import matplotlib.pyplot as plt
import os

dis_dir_unsorted = glob.glob('Ecoul_*')
list_dis_dirs = sorted(dis_dir_unsorted, key = lambda x: int(x.split('_')[-1]))
for dis in list_dis_dirs[:]:
    if os.path.exists(dis+"/analysis/act_energy/act_energy_DMR_set.txt"):
        act_energy = (np.loadtxt(dis+"/analysis/act_energy/act_energy_DMR_set.txt",comments='#',unpack=True))[1]
        dmr        = (np.loadtxt(dis+"/analysis/act_energy/act_energy_DMR_set.txt",comments='#',unpack=True))[0]
        dmr_dir_unsorted = glob.glob('DMR_*')
        list_DMR_dirs    = sorted(dmr_dir_unsorted, key = lambda x: int(x.split('_')[-1]))
        Ecoul         = (np.loadtxt(dis+"/analysis/sim_param_DMR_set.txt")[4])
        plt.plot(dmr*100,act_energy,marker="+", linestyle="None",label="$E_C$ = {} eV".format(Ecoul))
plt.xlabel('DMR (%)')
plt.xscale("log")
plt.ylabel('$E_A$ (meV)')
plt.legend()
plt.savefig('act_energy_CTenergyset.png')
plt.close()       