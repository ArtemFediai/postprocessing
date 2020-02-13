#!/usr/bin/env python3
import argparse
from SimAnalysis import CurrTempSimulation, CurrTempDMRset
import glob
import numpy as np
import matplotlib.pyplot as plt
import os

parser = argparse.ArgumentParser()

# parser.add_argument('-single',action='store_true', help='Specify if one single spec_temp set is analysed.')
parser.add_argument('-marcus', action='store_true', help='Specify if simulation was done with Marcus rates.')
parser.add_argument('-Ecoul', action='store_true', help='Specify if sim. was done for a specific Ecoul.')
parser.add_argument('-t','--spec_temp', type=int, help='Specify T for which to consider conductivity.')


args = parser.parse_args()
print(args)

if args.Ecoul:
    rates = "Marcus"  # currently implemented that way in SimAnalysis
else:
    if args.marcus:
        rates = "Marcus" 
    else:
        rates = "Miller"

spec_temp = args.spec_temp
print(spec_temp)

DMR_dir_unsorted = glob.glob('analysis/DMR_*')
list_DMR_dirs = sorted(DMR_dir_unsorted, key = lambda x: int(x.split('_')[-1]))

dmr_set = CurrTempDMRset(rates=rates, analyse_Ecoul=args.Ecoul)
# dmr_set.plot_av_current(Tlim_low=250, Tlim_high=400,errorbar=True)
dmr_set.plot_av_conductivity(errorbar=True,plot_log=True)

dmr          = dmr_set.DMR
n_DMR = len(list_DMR_dirs)
cond = np.zeros(n_DMR)
std_cond = np.zeros(n_DMR)
log_std_cond = np.zeros(n_DMR)



for i_dmr, dmr_dir in enumerate(list_DMR_dirs):
    cond_data = np.genfromtxt(dmr_dir+'/av_conductivity.txt', dtype = float)
    temp = cond_data[:,0]
    temp_idx = (np.abs(temp-spec_temp)).argmin()
    print("Found {} K. Use spec_temp = {} K from now on.".format(temp[temp_idx],temp[temp_idx]))    
    spec_temp = temp[temp_idx]
    cond[i_dmr]         = cond_data[temp_idx,1]
    std_cond[i_dmr]     = cond_data[temp_idx,2]
    log_std_cond[i_dmr] = cond_data[temp_idx,3]

# save data
if not os.path.isdir('analysis/cond_vs_dmr'):
    os.makedirs('analysis/cond_vs_dmr')
with open("analysis/cond_vs_dmr/cond_vs_dmr.txt",'w') as f:
    if rates == "Marcus":   
        if args.Ecoul:
            f.write('# Temp. = {} K, field = {} V/nm, disorder = {} eV, lambda = {} eV, Ecoul = {} eV\n'.format(spec_temp, dmr_set.field,dmr_set.dis,dmr_set.lam,dmr_set.Ecoul))
        else:
            f.write('# Temp. = {} K, field = {} V/nm, disorder = {} eV, lambda = {} eV\n'.format(spec_temp, dmr_set.field,dmr_set.dis,dmr_set.lam))
    else:    
        f.write('# Temp. = {} K, field = {} V/nm, disorder = {} eV\n'.format(spec_temp, dmr_set.field,dmr_set.dis))  
    f.write('#\n# DMR       Av. Conductivity(S/m)      Normal std. error(S/m)       Log. std. error(S/m)\n')
    for i in range(n_DMR):
        if cond[i]==0.0:
            continue
        else:
            f.write('  {0:<7}   {1:<24}   {2:<26}   {3:<1}\n'.format(dmr[i],cond[i],std_cond[i],log_std_cond[i]))
f.close()


# plot 
dmr = 100.0*np.array(dmr)
plot_log = True
errorbar = False


if plot_log:
    if errorbar:
        plt.errorbar(dmr,np.log(cond),yerr=log_std_cond)
    else:
        plt.plot(dmr,np.log(cond),marker="+", linestyle="None")
    ylabel = 'log $\sigma$ (S/m)'
else:
    if errorbar:
        plt.errorbar(dmr,cond,yerr=std_cond)
    else:
        plt.plot(dmr,cond,marker="+", linestyle="None")
    ylabel = '$\sigma$ (S/m)'
plt.xlabel('DMR ($\%$)')
plt.xscale('log')
plt.ylabel(ylabel)
# plt.legend()
if rates == "Marcus":   
    if args.Ecoul:
        plt.title('Temp. = {} K, field = {} V/nm, disorder = {} eV, lambda = {} eV, Ecoul = {} eV\n'.format(spec_temp, dmr_set.field,dmr_set.dis,dmr_set.lam,dmr_set.Ecoul))
    else:
        plt.title('Temp. = {} K, field = {} V/nm, disorder = {} eV, lambda = {} eV\n'.format(spec_temp, dmr_set.field,dmr_set.dis,dmr_set.lam))
else:    
    plt.title('Temp. = {} K, field = {} V/nm, disorder = {} eV\n'.format(spec_temp, dmr_set.field,dmr_set.dis))  
plt.savefig('analysis/cond_vs_dmr/cond_vs_dmr.png')
plt.close()