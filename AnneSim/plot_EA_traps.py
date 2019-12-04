#!/usr/bin/env python3
import numpy as np 
import glob
import matplotlib.pyplot as plt
import os



# # without trap
# without_outputdir = "no_trap/act_energy_DMR_set_all.txt"
# act_energy = (np.loadtxt(without_outputdir,comments='#',unpack=True))[1]
# dmr        = (np.loadtxt(without_outputdir,comments='#',unpack=True))[0]
# plt.plot(dmr*100,act_energy,marker="+", linestyle="None",label="no traps")
# traps
list_dirs = sorted(glob.glob('trap_*'), key = lambda x: int(x.split('_')[-1]))
for ll in list_dirs[:]:
    outputdir = ll+"/analysis/act_energy/act_energy_DMR_set.txt"
    if os.path.exists(outputdir):
        act_energy = (np.loadtxt(outputdir,comments='#',unpack=True))[1]
        dmr        = (np.loadtxt(outputdir,comments='#',unpack=True))[0]
        dmr_dir_unsorted = glob.glob('DMR_*')
        list_DMR_dirs = sorted(dmr_dir_unsorted, key = lambda x: int(x.split('_')[-1]))
        if ll == "trap_0":
            field, disorder = (np.loadtxt(ll+'/analysis/sim_param_DMR_set.txt')[1:3])
        trap_conc       = np.loadtxt(ll+'/analysis/sim_param_DMR_set.txt')[4]
        plt.plot(dmr*100,act_energy,marker="+", linestyle="None",label="trap conc. = {} %".format(100*trap_conc))
plt.title('Rates = {}, Field = {} eV, disorder = {} eV.'.format("mixed-marcus", field,disorder))
plt.xlabel('DMR (%)')
plt.xscale("log")
plt.ylabel('$E_A$ (meV)')
plt.legend()
plt.savefig('act_energy_traps_compare.png')
plt.close()       