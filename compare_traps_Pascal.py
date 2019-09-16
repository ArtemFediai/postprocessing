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

Names = ['TCTA', 'aNPD',  'TPBI']

for i in range(len(Names)):
    mobility = np.loadtxt('{}/m/postprocessing/mobility.dat'.format(Names[i]))  # in CI
    n = np.loadtxt('{}/m/postprocessing/n_traps.dat'.format(Names[i]))  # trap molar ratio
    plt.plot(n[0:len(mobility)], mobility*1E4, label=Names[i], marker='o', LineStyle='None')
plt.xscale('log')
plt.yscale('log')
plt.ylabel('mobility, cm$^2$V$^{-1}$s$^{-1}$')
plt.xlabel('trap molar fraction')
plt.legend()
plt.xlim([1E-6, 2E-4])
plt.ylim([1E-7, 3E-3])
axes = plt.gca()
ylims = axes.get_ylim()
plt.plot([1E-4, 1E-4], [ylims[0], ylims[1]], '--')
plt.savefig("mobility_{}.png".format(Names[:]), dpi=600)
plt.close()