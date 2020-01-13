#!/usr/bin/env python3
import numpy as np 
import glob
import matplotlib.pyplot as plt
import os
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-dis', action='store_true', help='Specify if EA vs. disorder should be determined.')
parser.add_argument('-lam', action='store_true', help='Specify if EA vs. lambda should be determined.')
args = parser.parse_args()
print(args)

if args.dis:
    dep = "dis"
    xlabel = "$\sigma$ (eV)"
elif args.lam:
    dep = "l"
    xlabel = "$\lambda$ (eV)"

list_dirs = ['cubic', 'real_morph']
for dirr in list_dirs:
    act_energy = (np.loadtxt(dirr+"/analysis/Eact_vs_{}.txt".format(dep),comments='#',unpack=True))[1]
    x          = (np.loadtxt(dirr+"/analysis/Eact_vs_{}.txt".format(dep),comments='#',unpack=True))[0]
    plt.plot(x,act_energy,marker="+", linestyle="None",label=dirr)

        
plt.xlabel(xlabel)
plt.ylabel('$E_A$ (meV)')
plt.legend()
plt.savefig('act_energy_{}_compare_real.png'.format(dep))
plt.close()       