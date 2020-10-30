import numpy as np
import matplotlib
from numpy.core._multiarray_umath import ndarray
import warnings
import inspect
import pickle

#matplotlib.use('pdf')
import tkinter
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import os
import re
import yaml
import scipy.constants as const
import time
import gzip

from IP_or_EA_vs_R_ver_4 import add_inverse_axis

##
def eps(radii, energies, ab,a1b1):
    c = 7.1997176734999995  # coefficient in the formula

    a, b, a1, b1 = ab[0], ab[1], a1b1[0], a1b1[1]  # analyze interval; plot interval

    target_i = np.where(np.logical_and(radii > a, radii < b))[0]
    # target_i1 = np.where(np.logical_and(radii > a1, radii < b1))[0]
    selected_r = radii[target_i]  # selected interval to compute epsilon
    # selected_r1 = radii[target_i1]
    selected_energies = energies[target_i]
    # selected_energies1 = energies[target_i1]

    coef_poly = np.polyfit(1.0 / selected_r, selected_energies, 1)

    selected_fitted_energies = coef_poly[0] / selected_r + coef_poly[1]
    # selected_fitted_energies1 = coef_poly[0] / selected_r1 + coef_poly[1]

    epsilon = c / (c + coef_poly[0])

    print( "eps = ", epsilon)

    e_a1 = coef_poly[0]/a1 + coef_poly[1]

    point_a = [0, coef_poly[1]]
    point_b = [10 / selected_r[-1], selected_energies[-1]]
    point_b1 = [10 / a1, e_a1]

    return point_a, point_b, point_b1, target_i
##


plt.figure('polarization average')

folders = ['C60' ,'aNPD' , 'TCTA']

radii = np.loadtxt(folders[0] + '/Analysis/p_elementwise/radii.dat')

for i, folder in enumerate(folders):
    e0i_anion = np.loadtxt(folder + '/Analysis/p_elementwise/dEe_anion.dat')  #DFT depol
    eii_anion = np.loadtxt(folder + '/Analysis/p_elementwise/dV0_anion.dat')  #class pol
    v0i_anion = np.loadtxt(folder + '/Analysis/p_elementwise/dVe_anion.dat')  #class depol
    e0i_cation = np.loadtxt(folder + '/Analysis/p_elementwise/dEe_cation.dat')  #DFT depol
    eii_cation = np.loadtxt(folder + '/Analysis/p_elementwise/dV0_cation.dat')  #class pol
    v0i_cation = np.loadtxt(folder + '/Analysis/p_elementwise/dVe_cation.dat')  #class depol

    p_av = (e0i_anion + eii_anion + v0i_anion + e0i_cation + eii_cation + v0i_cation)/2.0

    point_a, point_b, point_b1, target_i = eps(radii, p_av, [19, 41], [50, 1E100])

    plt.plot([point_a[0], point_b1[0]],[point_a[1],point_b1[1]], color='black')

    plt.plot(10*radii[target_i]**(-1), p_av[target_i], 'o', color='black', mfc='none')

    plt.plot(10*radii**(-1), p_av, label = folder)

plt.xlim([0, 2])
plt.grid()
plt.ylabel('Average polarization energy, eV')
plt.xlabel('$10/R,  10 / \AA^{-1} $')
plt.xlabel('Inverse distance $1/R,  10 / \AA^{-1} $')
plt.legend()
ax1 = plt.gca()
add_inverse_axis(ax1)


plt.savefig('p_average_for_all.png', dpi=600)

