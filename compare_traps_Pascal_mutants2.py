"""
To be used in the folder where you compute something.
Compare many-body and stupid DOS
"""
from __future__ import print_function
from __future__ import absolute_import
import extract
import numpy as np
import matplotlib
matplotlib.use('pdf')
import matplotlib.pyplot as plt
import os
from shutil import copyfile

logscale = False

Names = ['M6_','M5_']

n_r = 30

for i in range(len(Names)):
    mobility = np.loadtxt('{}/postprocessing/mobility.dat'.format(Names[i]))  # in CI
    std_err_mo = np.loadtxt('{}/postprocessing/std_err_mo.dat'.format(Names[i]))  # in CI


    n = np.loadtxt('{}/postprocessing/n_traps.dat'.format(Names[i]))  # trap molar ratio

    print("size mobility", len(mobility))
    print("size std mobility", len(std_err_mo))
    print("size n", len(n))


    #plt.errorbar(n[0:len(mobility)], mobility*1E4, yerr=std_err_mo*1E4*np.sqrt(n_r), label=Names[i], capsize=3, alpha=1, Color='C{}'.format(i), fmt='-', ecolor='C{}'.format(i))
    plt.errorbar(n[0:len(mobility)], mobility*1E4, yerr=std_err_mo*1E4*np.sqrt(n_r), label=Names[i], capsize=3, alpha=1, Color='C{}'.format(i), fmt='-', ecolor='grey')
    #plt.plot(n[0:len(mobility)], mobility*1E4, label=Names[i], marker='o', LineStyle='None', Color="C{}".format(i), mfc="C{}".format(i), mec=None)


Names = ['TCTA', 'aNPD',  'TPBI']

for i in range(len(Names)):
    mobility = np.loadtxt('/home/artem/Desktop/fh2_ws/2019/092019/04092019/{}/m/postprocessing/mobility.dat'.format(Names[i]))  # in CI
    std_err_mo = np.loadtxt('/home/artem/Desktop/fh2_ws/2019/092019/04092019/{}/m/postprocessing/std_err_mo.dat'.format(Names[i]))  # in CI
    n = np.loadtxt('/home/artem/Desktop/fh2_ws/2019/092019/04092019/{}/m/postprocessing/n_traps.dat'.format(Names[i]))  # trap molar ratio
    #plt.errorbar(n[0:len(mobility)], mobility*1E4, yerr=std_err_mo*1E4*np.sqrt(n_r), label=Names[i], capsize=3, alpha=1, Color='C{}'.format(i), fmt='-', ecolor='C{}'.format(i))
    plt.errorbar(n[0:len(mobility)], mobility*1E4, yerr=std_err_mo*1E4*np.sqrt(n_r), label=Names[i], capsize=3, alpha=1, Color='C{}'.format(i), fmt='-', ecolor='grey', LineStyle=':')
    #plt.plot(n[0:len(mobility)], mobility*1E4, label=Names[i], marker='o', LineStyle='None', Color="C{}".format(i), mfc="C{}".format(i), mec=None)



plt.xscale('log')
plt.yscale('log')
plt.ylabel('Electron mobility, cm$^2$V$^{-1}$s$^{-1}$')
plt.xlabel('Trap molar fraction')
plt.legend()
plt.xlim([1E-6, 2E-4])
plt.ylim([1E-7, 3E-3])
axes = plt.gca()
ylims = axes.get_ylim()
plt.plot([1E-4, 1E-4], [ylims[0], ylims[1]], '--', Color='C{}'.format(3))
plt.savefig("mobility_{}_grey.png".format(Names[:]), dpi=600)
plt.close()