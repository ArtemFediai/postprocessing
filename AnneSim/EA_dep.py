#!/usr/bin/env python3
import argparse
from SimAnalysis import CurrTempSimulation, CurrTempDMRset
import glob
import numpy as np
import matplotlib.pyplot as plt
import os



parser = argparse.ArgumentParser()
parser.add_argument('-ect', action='store_true', help='Specify if EA vs. Ect should be determined.')
parser.add_argument('-dis', action='store_true', help='Specify if EA vs. disorder should be determined.')
parser.add_argument('-lam', action='store_true', help='Specify if EA vs. lambda should be determined.')
args = parser.parse_args()
print(args)

rates = "Marcus" 

if args.ect:
    dep = "Ecoul"
elif args.dis:
    dep = "dis"
elif args.lam:
    dep = "l"

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
plt.title(sim_info,fontsize=14)
plt.savefig("analysis/Eact_vs_{}.png".format(dep))
plt.close()        
