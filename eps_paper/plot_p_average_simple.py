import numpy as np
import matplotlib
from numpy.core._multiarray_umath import ndarray
import warnings
import inspect
import pickle
matplotlib.use('pdf')
import matplotlib.pyplot as plt
import os
import re
import yaml
import scipy.constants as const
import time
import gzip

from IP_or_EA_vs_R_ver_4 import add_inverse_axis

##

def main():


    # 1. place this script outside the folder with data
    # 2. specify the folder with a yaml file
    # 3. specify a yaml file
    # enjoy

    folder = 'C60'
    filename = 'jonas_dict.yaml'
    filename = 'C60.yaml'

    with open(folder + '/' + filename) as fid:
        jonas_dict = yaml.load(fid, Loader=yaml.FullLoader)

    radii = np.array(jonas_dict['R_inv_grid'])**(-1)


    e0i = np.array(jonas_dict['E0_dft'])
    v0i = np.array(jonas_dict['E0_sd'])
    eii = np.array(jonas_dict['E0_env_env'])

    e0i_anion = e0i - np.array(jonas_dict['E_anion_dft'])
    v0i_anion = v0i - np.array(jonas_dict['E_anion_sd'])
    eii_anion = eii - np.array(jonas_dict['E_anion_env_env'])

    e0i_cation = -np.array(jonas_dict['E_cation_dft']) + e0i
    v0i_cation = -np.array(jonas_dict['E_cation_sd']) + v0i
    eii_cation = -np.array(jonas_dict['E_cation_env_env']) + eii

    p_av = (e0i_anion + eii_anion + v0i_anion + e0i_cation + eii_cation + v0i_cation) / 2.0

    eps_only_computed = eps_only(folder, radii, p_av, [20.0, 40.0])
    print("Computed eps with eps_only for " + folder, eps_only_computed)
    exit()




########################################################################################################################

    point_a, point_b, point_b1, target_i = eps(folder, radii, p_av, [20, 40], [50, 1E100])

    # for average eps -->
    plt.figure('polarization average')
    plt.xlim([0, 2])
    plt.grid()
    plt.ylabel('Average polarization energy, eV')
    plt.xlabel('$10/R,  10 / \AA^{-1} $')
    plt.xlabel('Inverse distance $1/R,  10 / \AA^{-1} $')
    plt.legend()
    ax1 = plt.gca()
    add_inverse_axis(ax1)
    plt.savefig('p_average_' + folder, dpi=600)


########################################### FUNCTIONS ##################################################################

def eps_only(folder, radii, energies, ab):

    c = 7.1997176734999995  # coefficient in the formula
    a, b = ab[0], ab[1]  # analyze interval; plot interval

    target_i = np.where(np.logical_and(radii >= a, radii <= b))[0]
    selected_r = radii[target_i]
    selected_energies = energies[target_i]
    coef_poly = np.polyfit(1.0 / selected_r, selected_energies, 1)
    epsilon = c / (c + coef_poly[0])
    return epsilon

def eps(folder, radii, energies, ab,a1b1):
    c = 7.1997176734999995  # coefficient in the formula

    a, b, a1, b1 = ab[0], ab[1], a1b1[0], a1b1[1]  # analyze interval; plot interval

    target_i = np.where(np.logical_and(radii >= a, radii <= b))[0]
    # target_i1 = np.where(np.logical_and(radii > a1, radii < b1))[0]
    selected_r = radii[target_i]  # selected interval to compute epsilon
    # selected_r1 = radii[target_i1]
    selected_energies = energies[target_i]
    # selected_energies1 = energies[target_i1]

    coef_poly = np.polyfit(1.0 / selected_r, selected_energies, 1)

    selected_fitted_energies = coef_poly[0] / selected_r + coef_poly[1]

    epsilon = c / (c + coef_poly[0])

    print( "eps " + folder + " = ", epsilon)

    e_a1 = coef_poly[0]/a1 + coef_poly[1]

    point_a = [0, coef_poly[1]]
    point_b = [10 / selected_r[-1], selected_energies[-1]]
    point_b1 = [10 / a1, e_a1]

    return point_a, point_b, point_b1, target_i

##
if __name__ == '__main__':
    main()