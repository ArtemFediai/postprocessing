from __future__ import print_function
from __future__ import absolute_import
import numpy as np
import matplotlib
matplotlib.use('pdf')
import matplotlib.pyplot as plt
import extract
import os

# extract Ion. dop. frac.
n_dis     = int(np.loadtxt('s4s.txt')[0])
n_temp    = int(np.loadtxt('s4s.txt')[1])
n_replica = int(np.loadtxt('s4s.txt')[2])
n_doMR    = 1
sim_data=np.genfromtxt('sim_param.txt',comments='#',dtype=float)
# field=sim_data[0]
temperature = sim_data[1:1+n_temp]
disorder    = sim_data[1+n_temp:1+n_temp+n_dis]
print(disorder)
# doping_ratio=sim_data[1+n_temp+n_dis:1+n_temp+n_dis+n_doMR]
count_nan = 0
if count_nan:
    savename = 'ion_dop_frac_fig5a_withnan.txt'
else:
    savename = 'ion_dop_frac_fig5a.txt'
if not os.path.exists(savename):   
    ion_dop_frac_av   = np.zeros((n_dis,n_temp))
    std_ion_dop_frac  = np.zeros((n_dis,n_temp))
    for i_dis in range(n_dis):
        for i_temp in range(n_temp):
            ion_dop_frac_r = []
            num_nan=0
            for i_r in range(n_replica):
                fn = 'jobs/dis_{}/temp_{}/r_{}/output_job_0'.format(i_dis,i_temp,i_r)
                print(fn)
                extract_ion_dop_frac = extract.extract_float(fn, "Ionized dopant fraction: ")
                if type(extract_ion_dop_frac)==float:
                    if np.log(extract_ion_dop_frac)>=0:
                        ion_dop_frac_r.append(extract.extract_e(fn, "Ionized dopant fraction: "))
                    else:
                        ion_dop_frac_r.append(extract_ion_dop_frac)
                else:
                    num_nan+=1
                    if count_nan:
                        ion_dop_frac_r.append(0.0)
                    
            print(ion_dop_frac_r)
            print('{} times nan encountered in {}'.format(num_nan,fn))
            ion_dop_frac_av[i_dis,i_temp]   = np.mean(ion_dop_frac_r)
            std_ion_dop_frac[i_dis,i_temp]  = np.std(ion_dop_frac_r)/np.sqrt(len(ion_dop_frac_r))
            print("Av. ionized dopant fraction: {}".format(ion_dop_frac_av[i_dis,i_temp]))
    # save DOS data
    with open(savename,'w') as f:
        f.write('# Ionized dopant fraction for disorder ')
        for i_dis in range(n_dis):
            f.write('{}, '.format(disorder[i_dis]))
        
        f.write('\n')   
        for i_dis in range(n_dis):
            for i_temp in range(n_temp):
                f.write('{} '.format(ion_dop_frac_av[i_dis, i_temp]))
            f.write('\n')
        f.write('# Standarderror\n')
        for i_dis in range(n_dis):
            for i_temp in range(n_temp):
                f.write('{} '.format(std_ion_dop_frac[i_dis, i_temp]))
            f.write('\n')
    f.close()




# plor Ion. dop. frac. vs. Offset and get simulation points
ion_dop_frac     = np.loadtxt(savename)[0:n_dis]
std_ion_dop_frac = np.loadtxt(savename)[n_dis:]
print(ion_dop_frac)

# matplotlib.rcParams.update({'font.size': 16})
# plt.figure(figsize = (5,5))
# plt.subplots_adjust(left=0.2,bottom=0.15)
for i_dis in range(n_dis):
    # plt.errorbar(1000/temperature,     ion_dop_frac[i_dis,:], yerr=std_ion_dop_frac[i_dis,:], lw=2, label='disorder= {} eV'.format(disorder[i_dis])) 
    if i_dis==0:
        lim=10
        plt.plot(1000/temperature[lim:],     ion_dop_frac[i_dis,lim:], lw=2, linestyle='None', marker='o', label='disorder= {} eV'.format(disorder[i_dis]))
    else:
        plt.plot(1000/temperature,     ion_dop_frac[i_dis,:], lw=2, linestyle='None', marker='o',label='disorder= {} eV'.format(disorder[i_dis])) 
plt.ylabel('Ionized dopant fraction')
plt.xlabel('Temperature [1000/K]')
plt.ylim(1e-04,1)
plt.yscale('log')
plt.legend()
if count_nan:
    plt.savefig('fig5a_withnan.png')
    plt.savefig('fig5a_withnan.svg')
else:
    plt.savefig('fig5a.png')
    plt.savefig('fig5a.svg')
plt.clf()
