#!/usr/bin/env python3
import numpy as np 
import glob
import matplotlib.pyplot as plt
import os

list_dirs = glob.glob('with*')
for ll in list_dirs[:]:
    outputdir = ll+"/act_energy_DMR_set_all.txt"
    if os.path.exists(outputdir):
        act_energy = (np.loadtxt(outputdir,comments='#',unpack=True))[1]
        dmr        = (np.loadtxt(outputdir,comments='#',unpack=True))[0]
        dmr_dir_unsorted = glob.glob('DMR_*')
        list_DMR_dirs = sorted(dmr_dir_unsorted, key = lambda x: int(x.split('_')[-1]))
        if ll == "with":
            field, disorder = (np.loadtxt(ll+'/low_DMR/analysis/sim_param_DMR_set.txt')[1:3])
        plt.plot(dmr*100,act_energy,marker="+", linestyle="None",label=ll)
plt.title('Rates = {}, Field = {} eV, disorder = {} eV.'.format("Marcus", field,disorder))
plt.xlabel('DMR (%)')
plt.xscale("log")
plt.ylabel('$E_A$ (meV)')
plt.legend()
plt.savefig('act_energy_traps_compare.png')
plt.close()       