#!/usr/bin/env python3
import numpy as np 
import glob
import matplotlib.pyplot as plt
import os

l_dir_unsorted = glob.glob('l_*')
list_l_dirs = sorted(l_dir_unsorted, key = lambda x: int(x.split('_')[-1]))
for ll in list_l_dirs[:]:
    if os.path.exists(ll+"/analysis/act_energy/act_energy_DMR_set.txt"):
        act_energy = (np.loadtxt(ll+"/analysis/act_energy/act_energy_DMR_set.txt",comments='#',unpack=True))[1]
        dmr        = (np.loadtxt(ll+"/analysis/act_energy/act_energy_DMR_set.txt",comments='#',unpack=True))[0]
        dmr_dir_unsorted = glob.glob('DMR_*')
        list_DMR_dirs = sorted(dmr_dir_unsorted, key = lambda x: int(x.split('_')[-1]))
        lam        = (np.loadtxt(ll+'/analysis/sim_param_DMR_set.txt')[3])
        if ll == "l_4":
            field, disorder = (np.loadtxt(ll+'/analysis/sim_param_DMR_set.txt')[1:3])
        plt.plot(dmr*100,act_energy,marker="+", linestyle="None",label="$\lambda$ = {} eV".format(lam))
plt.title('Rates = {}, Field = {} eV, disorder = {} eV.'.format("Marcus", field,disorder))
plt.xlabel('DMR (%)')
plt.xscale("log")
plt.ylabel('$E_A$ (meV)')
plt.legend()
plt.savefig('act_energy_lambdaset.png')
plt.close()       