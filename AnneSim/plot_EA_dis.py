#!/usr/bin/env python3
import numpy as np 
import glob
import matplotlib.pyplot as plt
import os

dis_dir_unsorted = glob.glob('dis_*')
list_dis_dirs = sorted(dis_dir_unsorted, key = lambda x: int(x.split('_')[-1]))
num_sim=2
for dis in list_dis_dirs[:]:
    if os.path.exists(dis+"/"+str(num_sim)+"/analysis/act_energy/act_energy_DMR_set.txt"):
        act_energy = (np.loadtxt(dis+"/"+str(num_sim)+"/analysis/act_energy/act_energy_DMR_set.txt",comments='#',unpack=True))[1]
        dmr        = (np.loadtxt(dis+"/"+str(num_sim)+"/analysis/act_energy/act_energy_DMR_set.txt",comments='#',unpack=True))[0]
        dmr_dir_unsorted = glob.glob('DMR_*')
        list_DMR_dirs    = sorted(dmr_dir_unsorted, key = lambda x: int(x.split('_')[-1]))
        disorder         = (np.loadtxt(dis+"/"+str(num_sim)+'/analysis/sim_param_DMR_set.txt')[2])
        plt.plot(dmr*100,act_energy,marker="+", linestyle="None",label="$\sigma$ = {} eV".format(disorder))
plt.xlabel('DMR (%)')
plt.xscale("log")
plt.ylabel('$E_A$ (meV)')
plt.legend()
plt.savefig('act_energy_disorderset_{}.png'.format(num_sim))
plt.close()       