"""
To be used in the folder where you compute something.
First plot: mobility vs. trap molar fraction. Red line: molar fraction = 1E-4
Second plot: mobility vs. [traps per electron]
"""
from __future__ import print_function
from __future__ import absolute_import
import numpy as np
import matplotlib
matplotlib.use('pdf')
import matplotlib.pyplot as plt
import os
from shutil import copyfile

logscale = False
Names = ['aNPD',  'TCTA', 'TPBI']
n_r = 30  # num of replicas


# First plot: mobility vs trap molar fraction

for i in range(len(Names)):
    mobility = np.loadtxt('{}/m/postprocessing/mobility.dat'.format(Names[i]))  # in CI
    std_err_mo = np.loadtxt('{}/m/postprocessing/std_err_mo.dat'.format(Names[i]))  # in CI
    n = np.loadtxt('{}/m/postprocessing/n_traps.dat'.format(Names[i]))  # trap molar ratio
    #plt.errorbar(n[0:len(mobility)], mobility*1E4, yerr=std_err_mo*1E4*np.sqrt(n_r), label=Names[i], capsize=3, alpha=1, Color='C{}'.format(i), fmt='-', ecolor='C{}'.format(i))
    plt.errorbar(n[0:len(mobility)], mobility*1E4, yerr=std_err_mo*1E4*np.sqrt(n_r), label=Names[i], capsize=3, alpha=1, Color='C{}'.format(i), fmt='-', ecolor='grey')
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


# Second plot: mobility vs [traps per electron]

density = np.array([1.131751335,	0.8779955674,	1.039569444])  # Deposit density [molecules*nm**-3]
num_sites = density*1E6  # 100x100x100
n_el = 100  # number of electrons
coef = num_sites/n_el  # number of sites per electron, for each mater.

for i in range(len(Names)):
    mobility = np.loadtxt('{}/m/postprocessing/mobility.dat'.format(Names[i]))  # in CI
    std_err_mo = np.loadtxt('{}/m/postprocessing/std_err_mo.dat'.format(Names[i]))  # in CI
    n = np.loadtxt('{}/m/postprocessing/n_traps.dat'.format(Names[i]))  # trap molar ratio as it is set
    n_t_over_n_el = coef[i]*n
    #plt.errorbar(n_t_over_n_el[0:len(mobility)], mobility*1E4, yerr=std_err_mo*1E4*np.sqrt(n_r), label=Names[i], capsize=3, alpha=1, Color='C{}'.format(i), fmt='-', ecolor='C{}'.format(i))
    plt.errorbar(n_t_over_n_el[0:len(mobility)], mobility*1E4, yerr=std_err_mo*1E4*np.sqrt(n_r), label=Names[i], capsize=3, alpha=1, Color='C{}'.format(i), fmt='-', ecolor='grey')
    #plt.plot(n[0:len(mobility)], mobility*1E4, label=Names[i], marker='o', LineStyle='None', Color="C{}".format(i), mfc="C{}".format(i), mec=None)

plt.xscale('log')
plt.yscale('log')
plt.ylabel('Electron mobility, cm$^2$V$^{-1}$s$^{-1}$')
plt.xlabel('Trap number over electron number')
plt.legend()
plt.xlim([2E-2, 2E0])
plt.ylim([1E-7, 3E-3])
axes = plt.gca()
ylims = axes.get_ylim()
plt.plot([1E-0, 1E-0], [ylims[0], ylims[1]], '--', Color='C{}'.format(3))
plt.savefig("mobility_{}_grey_NEW.png".format(Names[:]), dpi=600)
plt.close()