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
    radii_not_renormalized = np.loadtxt(folders[0] + '/Analysis/p_elementwise/radii.dat')


    # for average eps -->
    a, b, c, n = 20, 30, 10, 50
    lower_limits_eps = np.linspace(a, a+c, n)
    upper_limits_eps = np.linspace(b, b+c, n)

    # <-- for average eps
    number_density_nm3_aNPD = 1.158623303810171  # aNPD
    number_density_nm3_TCTA = 0.9162165860886232  # TCTA
    number_density_nm3_C60 = 1.4266923  # C60

    plt.figure('R vs R')

    name_of_postprocessing_folder = 'postprocessing_output'
    name_of_pcs_pcs_folder = 'pcs_pcs'
    name_of_pcs_qp_folder = 'pcs_qp'

    if not os.path.exists(name_of_postprocessing_folder):
        os.mkdir(name_of_postprocessing_folder)
        if not os.path.exists(name_of_postprocessing_folder + '/' + name_of_pcs_qp_folder):
            os.mkdir(name_of_postprocessing_folder + '/' + name_of_pcs_qp_folder)
    else:
        pass

    ab = [0, 2]
    plt.plot(ab, ab, color='grey', linestyle=':', label='homogenious material')
    for i, folder in enumerate(folders):
        radii = np.loadtxt(folder + '/Analysis/p_elementwise/radii_renormalized.dat')
        e0i_anion = np.loadtxt(folder + '/Analysis/p_elementwise/dEe_anion.dat')  #DFT depol
        eii_anion = np.loadtxt(folder + '/Analysis/p_elementwise/dV0_anion.dat')  #class pol
        v0i_anion = np.loadtxt(folder + '/Analysis/p_elementwise/dVe_anion.dat')  #class depol
        e0i_cation = np.loadtxt(folder + '/Analysis/p_elementwise/dEe_cation.dat')  #DFT depol
        eii_cation = np.loadtxt(folder + '/Analysis/p_elementwise/dV0_cation.dat')  #class pol
        v0i_cation = np.loadtxt(folder + '/Analysis/p_elementwise/dVe_cation.dat')  #class depol

        p_av = (e0i_anion + eii_anion + v0i_anion + e0i_cation + eii_cation + v0i_cation)/2.0

        point_a, point_b, point_b1, target_i = eps(folder, radii, p_av, [20, 40], [50, 1E100])

        #plt.plot([point_a[0], point_b1[0]], [point_a[1],point_b1[1]], color='black')
        #plt.plot(10*radii[target_i]**(-1), p_av[target_i], 'o', color='black', mfc='none')
        #plt.plot(10*radii**(-1), p_av, label=folder)


        plt.plot(10*radii_not_renormalized**-1, 10*radii**-1, label=folders_names[i], marker='o', markerfacecolor='None')

        #save in postprocessing folder
        #np.savetxt(name_of_postprocessing_folder + '/' + name_of_pcs_qp_folder + '/' + 'radii.dat', radii_not_renormalized)
        #np.savetxt(name_of_postprocessing_folder + '/' + name_of_pcs_qp_folder + '/' + 'p_av_{}.dat'.format(folder), p_av)

        #eps_range = np.zeros(n)
        #for i in range(n):
        #    eps_range[i] = eps_only(folder, radii, p_av, [lower_limits_eps[i], upper_limits_eps[i]])
        #eps_std = np.std(eps_range)
        #eps_mean = np.mean(eps_range)

        #print('eps_mean ' + folder + '= ', eps_mean)
        #print('eps_std ' + folder + '= ', eps_std)

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


    plt.draw()
    plt.close()

    #
    #ax, smth = plt.subplots()
    #ax.plot(radii, radii_not_renormalized)

    #



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

if __name__ == '__main__':
    main()