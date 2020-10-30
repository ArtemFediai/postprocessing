"""
IP(R) plot from different files
"""

from __future__ import print_function
from __future__ import absolute_import
import numpy as np
import matplotlib
matplotlib.use('pdf')
import matplotlib.pyplot as plt
import os
import re
import yaml
import scipy.constants as const
import time

basis_path_lib = {'def2-SVP': 'aNPD_def2_SVP', 'def2-SVPD': 'aNPD_basis_set', 'nm-SVPD': 'nm_SVPD'}
basis_IP_lib = {'def2-SVP': [], 'def2-SVPD': [], 'nm-SVPD': []}
basis_R_lib = {'def2-SVP': [], 'def2-SVPD': [], 'nm-SVPD': []}


len(basis_path_lib)

for key in basis_path_lib:
    IP = np.loadtxt(basis_path_lib[key] + "/IP.dat")
    R = np.loadtxt(basis_path_lib[key] + "/r.dat")
    basis_IP_lib[key] = IP
    basis_R_lib[key] = R

# IP vs R
for key in basis_path_lib:
    plt.plot(basis_R_lib[key], basis_IP_lib[key], LineStyle='-', marker='o')
plt.xlabel('R, A')
plt.ylabel('IP, eV')
plt.legend(basis_path_lib.keys())
plt.savefig('IP_vs_R_sd.png')
plt.close()

# IP vs 1/R
for key in basis_path_lib:
    plt.plot(10/basis_R_lib[key], basis_IP_lib[key], LineStyle='-', marker='o')
plt.xlabel('10/R, A')
plt.ylabel('IP, eV')
plt.legend(basis_path_lib.keys())
plt.savefig('Inverse_IP_vs_R_sd.png')
plt.close()

print("I am done")