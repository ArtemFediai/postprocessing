#!/usr/bin/env python3
import numpy as np 
import glob
import matplotlib.pyplot as plt
import os


dep    = "DMR"
xlabel = "DMF ($\%$)"

list_dirs = ['cubic', 'real_morph']
for dirr in list_dirs:
    act_energy = (np.loadtxt(dirr+"/analysis/act_energy/act_energy_DMR_set.txt",comments='#',unpack=True))[1]
    x          = (np.loadtxt(dirr+"/analysis/act_energy/act_energy_DMR_set.txt",comments='#',unpack=True))[0]
    plt.plot(x*100,act_energy,marker="+", linestyle="None",label=dirr)

        
plt.xlabel(xlabel)
plt.xscale("log")
plt.ylabel('$E_A$ (meV)')
plt.legend()
plt.savefig('act_energy_{}_compare_real.png'.format(dep))
plt.close()       