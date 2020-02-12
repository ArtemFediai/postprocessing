from __future__ import print_function
from __future__ import absolute_import
import numpy as np
import matplotlib
matplotlib.use('pdf')
import matplotlib.pyplot as plt
import extract
import os

# extract Ion. dop. frac.
n_offset  = int(np.loadtxt('s4s.txt')[0])
n_temp    = int(np.loadtxt('s4s.txt')[1])
n_replica = int(np.loadtxt('s4s.txt')[2])
# n_offset    = 4
sim_data=np.genfromtxt('sim_param.txt',comments='#',dtype=float)
# field=sim_data[0]
temperature = sim_data[1:1+n_temp]
offset    = sim_data[1+n_temp:1+n_temp+n_offset]
disorder  = float(sim_data[1+n_temp+n_offset:1+n_temp+n_offset+1])
# doping_ratio=sim_data[1+n_temp+n_offset:1+n_temp+n_offset+n_doMR]

if not os.path.exists('ion_dop_frac_fig5b.txt'):   
    ion_dop_frac_av  = np.zeros((n_offset,n_temp))
    std_ion_dop_frac = np.zeros((n_offset,n_temp))
    for i_off in range(n_offset):
        for i_temp in range(n_temp):
            ion_dop_frac_r = []
            count_nan = 0
            for i_r in range(n_replica):
                fn = 'jobs/off_{}/temp_{}/r_{}/output_job_0'.format(i_off,i_temp,i_r)
                print(fn)
                extract_ion_dop_frac = extract.extract_float(fn, "Ionized dopant fraction: ")
                if type(extract_ion_dop_frac)==float:
                    if np.log(extract_ion_dop_frac)>=0:
                        ion_dop_frac_r.append(extract.extract_e(fn, "Ionized dopant fraction: "))
                    else:
                        ion_dop_frac_r.append(extract.extract_float(fn, "Ionized dopant fraction: "))
                else:
                    count_nan+=1
            print('{} times nan encountered in {}'.format(count_nan,fn))
            ion_dop_frac_av[i_off,i_temp]  = np.mean(ion_dop_frac_r)
            std_ion_dop_frac[i_off,i_temp] = np.std(ion_dop_frac_r)/np.sqrt(len(ion_dop_frac_r))
            print("Av. ionized dopant fraction: {}".format(ion_dop_frac_av[i_off,i_temp]))

    # save DOS data
    with open('ion_dop_frac_fig5b.txt','w') as f:
        f.write('# Ionized dopant fraction for disorder = {} and offset = '.format(disorder))
        for i_off in range(n_offset):
            f.write('{}, '.format(offset[i_off]))
        f.write('\n')
        for i_off in range(n_offset):
            for i_temp in range(n_temp):
                f.write('{} '.format(ion_dop_frac_av[i_off, i_temp]))
            f.write('\n')
        f.write('# Standarderror\n')
        for i_dis in range(n_offset):
            for i_temp in range(n_temp):
                f.write('{} '.format(std_ion_dop_frac[i_off, i_temp]))
            f.write('\n')
    f.close()




# plor Ion. dop. frac. vs. Offset and get simulation points
ion_dop_frac     = np.loadtxt('ion_dop_frac_fig5b.txt')[0:n_offset]
std_ion_dop_frac = np.loadtxt('ion_dop_frac_fig5b.txt')[n_offset:]

# matplotlib.rcParams.update({'font.size': 16})
# plt.figure(figsize = (5,5))
# plt.subplots_adjust(left=0.2,bottom=0.15)
colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728',
              '#9467bd', '#e377c2']
for i_off in range(n_offset):
    # plt.errorbar(1000/temperature,     ion_dop_frac[i_off,:], yerr=std_ion_dop_frac[i_off,:], lw=2, label='offset= {} eV'.format(offset[i_off])) 
    plt.plot(1000/temperature,     ion_dop_frac[i_off,:], lw=2,color=colors[i_off], linestyle = 'None', marker='o', label='offset = {} eV'.format(-offset[i_off])) 
plt.ylabel('Ionized dopant fraction')
plt.xlabel('Temperature [1000/K]')
# plt.yscale('log')
plt.ylim(0,1)
plt.legend()
plt.savefig('fig5b.png')
plt.savefig('fig5b.svg')
plt.clf()
