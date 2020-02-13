from __future__ import print_function
from __future__ import absolute_import
import numpy as np
import matplotlib
matplotlib.use('pdf')
import matplotlib.pyplot as plt
import extract
import os

# extract Ion. dop. frac.
if not os.path.exists('ion_dop_frac.txt'):
    n_offset  = int(np.loadtxt('s4s.txt')[0])
    n_replica = int(np.loadtxt('s4s.txt')[1])

    
    ion_dop_frac_av = np.zeros(n_offset)
    for i_off in range(n_offset):
        ion_dop_frac_r = np.zeros(n_replica)
        for i_r in range(0,n_replica):
            fn = 'jobs/off_{}/r_{}/output_job_0'.format(i_off,i_r)
            print(fn)
            ion_dop_frac_r[i_r] = extract.extract_float(fn, "Ionized dopant fraction: ")
            if np.log(ion_dop_frac_r[i_r]) > 0:
                ion_dop_frac_r[i_r] = extract.extract_e(fn, "Ionized dopant fraction: ")
        ion_dop_frac_av[i_off]		= np.mean(ion_dop_frac_r)
        print(ion_dop_frac_r)
        print("Av. ionized dopant fraction: {}".format(ion_dop_frac_av[i_off]))

    # save DOS data
    with open('ion_dop_frac.txt','w') as f:
        f.write('# Ionized dopant fraction\n')
        for i_off in range(n_offset):
            f.write('{}\n'.format(ion_dop_frac_av[i_off]))
    f.close()




# plor Ion. dop. frac. vs. Offset and get simulation points
offset        = np.loadtxt('offset.txt')
ion_dop_frac  = np.loadtxt('ion_dop_frac.txt')

idx_0ion   = 0
idx_50ion  = 10
idx_100ion = 20

matplotlib.rcParams.update({'font.size': 16})
plt.figure(figsize = (5,5))
plt.subplots_adjust(left=0.2,bottom=0.15)
plt.plot(offset,           ion_dop_frac[:], lw=2) # TODO: Ask Artem which row (0 or 1) to take ?
msize = 5
plt.plot(offset[idx_0ion],   ion_dop_frac[idx_0ion],   label='{0:.2%}'.format(ion_dop_frac[idx_0ion]),marker='o', markersize=msize,color='b')
plt.plot(offset[idx_50ion],  ion_dop_frac[idx_50ion],  label='{0:.2%}'.format(ion_dop_frac[idx_50ion]),marker='o', markersize=msize,color='k')
plt.plot(offset[idx_100ion], ion_dop_frac[idx_100ion], label='{0:.2%}'.format(ion_dop_frac[idx_100ion]),marker='o', markersize=msize,color='r')
plt.grid()
plt.ylim(0,1.05)
plt.ylabel('Ionized dopant fraction')
plt.xlabel('Offset [eV]')
plt.legend()
plt.savefig('fig4a.png')
plt.savefig('fig4a.svg')
plt.clf()

# write to txt file
with open('fig4_simpoints.txt','w') as f:
    f.write('# fig 4\n')
    f.write('# Offset   Ion. dopant fraction\n')
    f.write('{} {}\n'.format(-offset[idx_0ion],   ion_dop_frac[idx_0ion]))
    f.write('{} {}\n'.format(-offset[idx_50ion],  ion_dop_frac[idx_50ion]))
    f.write('{} {}\n'.format(-offset[idx_100ion], ion_dop_frac[idx_100ion]))
f.close()