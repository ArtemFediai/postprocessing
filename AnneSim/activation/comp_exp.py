#!/usr/bin/env python3
import numpy as np 
import glob
import matplotlib.pyplot as plt
import os

# make plot nicer
msize = 5
mstyle = "o"

# exp points
x_cr2 = 100*np.array([0.002, 0.0033, 0.0041, 0.0054, 0.011, 0.022, 0.045, 0.097, 0.15, 0.21, 0.35])
y_cr2 = np.array([145, 90, 110, 85, 87, 75, 65, 60, 57, 75, 90])

x_w2 = 100*np.array([0.004, 0.0079, 0.008, 0.016, 0.033, 0.069, 0.15, 0.25, 0.37])
y_w2 = np.array([160, 115, 112, 99, 70, 75, 70, 79, 95])

# sim points
act_energy = (np.loadtxt("analysis/act_energy/act_energy_DMR_set.txt",comments='#',unpack=True))[1]
dmr        = (np.loadtxt("analysis/act_energy/act_energy_DMR_set.txt",comments='#',unpack=True))[0]
plt.plot(dmr*100,act_energy, marker=mstyle, markersize=msize, linestyle="None",label="sim")
plt.plot(x_cr2,y_cr2,label="C$_{60}$:Cr$_2$(hpp)$_4$",linestyle="None",marker="+")
plt.plot(x_w2,y_w2,label="C$_{60}$:W$_2$(hpp)$_4$",linestyle="None",marker="+")
plt.xlabel('Doping [%]')
plt.xscale("log")
plt.ylabel('Activation energy [meV]')
plt.legend()
plt.savefig('act_energy_comp_exp.png')
plt.savefig('act_energy_comp_exp.svg')
plt.savefig('act_energy_comp_exp.pdf')
plt.close()       