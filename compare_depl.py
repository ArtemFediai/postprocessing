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


Names = ['.', 'replay_100nm', 'replay_longer']

n_r = 30

for i in range(len(Names)):
    mobile_fraction = np.loadtxt('{}/postprocessing/mobile_fraction.dat'.format(Names[i]))
    std_err_mo = np.loadtxt('{}/postprocessing/std_err_mo.dat'.format(Names[i]))  # in CI
    dop = np.loadtxt('{}/postprocessing/dop.dat'.format(Names[i]))  # dop molar fraction
    plt.errorbar(dop[0:len(mobile_fraction)], mobile_fraction, yerr=std_err_mo*np.sqrt(n_r), label=Names[i], capsize=3, alpha=1, Color='C{}'.format(i), fmt='-', ecolor='C{}'.format(i))
    #plt.errorbar(dop[0:len(mobile_fraction)], mobile_fraction, yerr=std_err_mo*np.sqrt(n_r), label=Names[i], capsize=3, alpha=1, Color='C{}'.format(i), fmt='-', ecolor='grey')

#Fig 1
plt.xscale('log')
plt.ylabel('Mobile Carrier Fraction')
plt.xlabel('Doping')
plt.legend()
plt.xlim([1E-3, 1E-1])
plt.ylim([0, 1.2])
plt.savefig("CompareMobFraction.png".format(Names[:]), dpi=600)
plt.close()

# Fig 2 || Normalized std
for i in range(len(Names)):
    dop = np.loadtxt('{}/postprocessing/dop.dat'.format(Names[i]))  # dop molar fraction
    mobile_fraction = np.loadtxt('{}/postprocessing/mobile_fraction.dat'.format(Names[i]))
    std_err_mo = np.loadtxt('{}/postprocessing/std_err_mo.dat'.format(Names[i]))  # in CI
    plt.plot(dop[0:len(mobile_fraction)], std_err_mo*np.sqrt(n_r)/mobile_fraction, label=Names[i])
plt.xscale('log')
plt.ylabel('std Mobile Carrier Fraction')
plt.xlabel('Doping')
plt.legend()
plt.xlim([1E-3, 1E-1])
plt.savefig("CompareStdMobFraction.png".format(Names[:]), dpi=600)
plt.close()