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
import matplotlib.pyplot as plt

from IP_or_EA_vs_R_ver_4 import add_inverse_axis

##

def main():

    folders = ['C60' ,'aNPD' , 'TCTA']
    folders_names = ['C60' ,'NPB' , 'TCTA']


    # for average eps -->
    a, b, c, n = 20, 30, 10, 50
    lower_limits_eps = np.linspace(a, a+c, n)
    upper_limits_eps = np.linspace(b, b+c, n)


    plt.figure('R vs R')

    name_of_postprocessing_folder = 'postprocessing_output'
    name_of_pcs_pcs_folder = 'pcs_pcs'
    name_of_pcs_qp_folder = 'pcs_qp'

    p_av_for_all = []
    radii_renormalized_all = []
    ab = [0, 2]
    plt.plot(ab, ab, color='grey', linestyle=':', label='homogenious material')
    for i, folder in enumerate(folders):
        with open(folder+'/jonas_dict.yaml') as fid:
            jonas_dict = yaml.load(fid, Loader=yaml.FullLoader)
        radii_not_renormalized = np.array(jonas_dict['R_direct_grid'])
        radii_renormalized = np.array(jonas_dict['r_renormalized'])

        e0i = np.array(jonas_dict['E0_dft'])
        v0i = np.array(jonas_dict['E0_sd'])
        eii = np.array(jonas_dict['E0_env_env'])

        e0i_anion = e0i - np.array(jonas_dict['E_anion_dft'])
        v0i_anion = v0i - np.array(jonas_dict['E_anion_sd'])
        eii_anion = eii - np.array(jonas_dict['E_anion_env_env'])

        e0i_cation = -np.array(jonas_dict['E_cation_dft']) + e0i
        v0i_cation = -np.array(jonas_dict['E_cation_sd']) + v0i
        eii_cation = -np.array(jonas_dict['E_cation_env_env']) + eii


        p_av = (e0i_anion + eii_anion + v0i_anion + e0i_cation + eii_cation + v0i_cation)/2.0
        p_av_for_all.append(p_av)
        radii_renormalized_all.append(radii_renormalized)
        plt.plot(10*radii_not_renormalized**(-1), 10*radii_renormalized**(-1), label=folders_names[i], marker='o', markerfacecolor='None', color='C{}'.format(i))


    plt.xlim([0, 2])
    plt.ylim([0, 2])
    plt.grid()
    plt.xlabel('$10/R, \mathrm{\AA}^{-1}$')
    plt.ylabel(r'10/ $\tilde{R},  \mathrm{\AA}^{-1} $')
    plt.legend()
    #plt.tight_layout()
    # ax1 = plt.gca()
    # add_inverse_axis(ax1)
    plt.gca().set_aspect('equal', adjustable='box')
    plt.savefig('R_vs_R.png', dpi=600)
    plt.savefig('R_vs_R.svg')
    plt.close()

    # pol_vs_R and pol_vs_R

    fig, axes = plt.subplots(1,2, figsize=[10,5])
    ax = axes[1]
    ax0 = axes[0]

    ax0.plot(ab, ab, color='grey', linestyle=':', label='homogenious material')
    for ii, folder_name in enumerate(folders_names):

        point_a, point_b, point_b1, target_i, point_a1, epsilon = eps_mod(folders[ii], radii_renormalized_all[ii], radii_not_renormalized, p_av_for_all[ii], [20, 40], [1, 1E5])
        ax.plot([point_a1[0], point_b1[0]], [point_a1[1], point_b1[1]], color='black')
        ax.plot(10*radii_renormalized_all[ii]**(-1.0), p_av_for_all[ii], linestyle='solid', label='$\mathrm{P}$' + '$(R)$'.format(folder_name) + ' for ' + folder_name, color='C{}'.format(ii))
        ax.plot(10*radii_not_renormalized**(-1.0), p_av_for_all[ii], linestyle='--', label= '$\mathrm{P}$' + r'$(\tilde{R})$' + ' for ' + folder_name, color='C{}'.format(ii), alpha=0.5, linewidth='3')

        ax0.plot(10*radii_not_renormalized**(-1), 10*radii_renormalized_all[ii]**(-1), label=folders_names[ii], markerfacecolor='None', color='C{}'.format(ii))


    #ax
    ax.set_xlabel(r'$10/R$, $10/ \tilde{R}, \mathrm{\AA}^{-1}$')
    ax.set_ylabel('Average polarization energy $\mathrm{P}$, eV')
    ax.set_xlim([0,1])
    ax.set_ylim([0.0,1.0])
    ax.legend()
    add_inverse_axis(ax)


    #ax0

    ax0.set_xlim([0, 2])
    ax0.set_ylim([0, 2])
    ax0.grid()
    ax0.set_xlabel('$10/R, \mathrm{\AA}^{-1}$')
    ax0.set_ylabel(r'10/ $\tilde{R},  \mathrm{\AA}^{-1} $')
    ax0.legend()
    #plt.tight_layout()
    # ax1 = plt.gca()
    # add_inverse_axis(ax1)
    ax0.set_aspect('equal', adjustable='box')

    # savefig
    fig.tight_layout()
    fig.savefig('P_vs_R_vs_R.png', dpi=600)
    fig.savefig('P_vs_R_vs_R.svg')






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

def eps_only(folder, radii, energies, ab):

    c = 7.1997176734999995  # coefficient in the formula
    a, b = ab[0], ab[1]  # analyze interval; plot interval

    target_i = np.where(np.logical_and(radii >= a, radii <= b))[0]
    selected_r = radii[target_i]
    selected_energies = energies[target_i]
    coef_poly = np.polyfit(1.0 / selected_r, selected_energies, 1)
    epsilon = c / (c + coef_poly[0])
    return epsilon

def eps_mod(folder, radii, r_not_renormalized, energies, ab,a1b1):
    # differs from eps by outputting another points
    c = 7.1997176734999995  # coefficient in the formula

    a, b, a1, b1 = ab[0], ab[1], a1b1[0], a1b1[1]  # analyze interval; plot interval

    target_i = np.where(np.logical_and(r_not_renormalized >= a, r_not_renormalized <= b))[0]
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
    e_b1 = coef_poly[0]/b1 + coef_poly[1]

    point_a = [10 / b1, coef_poly[1]]
    point_b = [10 / selected_r[-1], selected_energies[-1]]
    point_b1 = [10 / a1, e_a1]
    point_a1 = [10 / b1, e_b1]

    return point_a, point_b, point_b1, target_i, point_a1, epsilon
if __name__ == '__main__':
    main()