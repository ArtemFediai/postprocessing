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

    folders = ['aNPD' , 'C60' ,'TCTA']
    radii = np.loadtxt(folders[0] + '/Analysis/p_elementwise/radii.dat')

    radii_renormalization = False


    # for average eps -->
    a, b, c, n = 20, 30, 10, 50
    lower_limits_eps = np.linspace(a, a+c, n)
    upper_limits_eps = np.linspace(b, b+c, n)

    # <-- for average eps
    number_density_nm3_aNPD = 1.158623303810171  # aNPD
    number_density_nm3_TCTA = 0.9162165860886232  # TCTA
    number_density_nm3_C60 = 1.4266923  # C60




    for i, folder in enumerate(folders):
        plt.figure('polarization_componentwise_{}'.format(folder), figsize=[4.5, 3.5])
        with open(folder+'/jonas_dict.yaml') as fid:
            jonas_dict = yaml.load(fid, Loader=yaml.FullLoader)
        radii_not_renormalized = np.array(jonas_dict['R_direct_grid'])
        if radii_renormalization == True:
            radii = np.array(jonas_dict['r_renormalized'])
        else:
            radii = np.array(jonas_dict['R_direct_grid'])

        e0i = np.array(jonas_dict['E0_dft'])
        v0i = np.array(jonas_dict['E0_sd'])
        eii = np.array(jonas_dict['E0_env_env'])
        e0 = np.array(jonas_dict['E0_E_core'])

        e0i_anion = e0i - np.array(jonas_dict['E_anion_dft'])
        v0i_anion = v0i - np.array(jonas_dict['E_anion_sd'])
        eii_anion = eii - np.array(jonas_dict['E_anion_env_env'])
        e0_anion = e0 - np.array(jonas_dict['E_anion_E_core'])

        e0i_cation = -np.array(jonas_dict['E_cation_dft']) + e0i
        v0i_cation = -np.array(jonas_dict['E_cation_sd']) + v0i
        eii_cation = -np.array(jonas_dict['E_cation_env_env']) + eii
        e0_cation = -np.array(jonas_dict['E_cation_E_core']) + e0

        delta_delta_anion = np.array(jonas_dict['delta_delta_E0_anion_mean'])
        delta_delta_cation = np.array(jonas_dict['delta_delta_E0_cation_mean'])

        p_av = (e0i_anion + eii_anion + v0i_anion + e0i_cation + eii_cation + v0i_cation)/2.0

        plt.plot(10*radii**(-1), 0.5*(e0i_anion+e0i_cation), label='$\Delta E_{env}$') #
        plt.plot(10*radii**(-1), 0.5*(v0i_anion+v0i_cation), label='$\Delta V_{core}$') #
        plt.plot(10*radii**(-1), 0.5*(eii_anion+eii_cation), label='$\Delta V_{env}$') #
#        plt.plot(10*radii**(-1), 0.5*(e0_anion+e0_cation)*np.ones(len(eii_cation)), label='$\Delta E_0$') #

        point_a, point_b, point_b1, target_i, point_a1, epsilon = eps_mod(folder, radii, radii_not_renormalized, p_av, [20, 40], [1, 1E5])


        # style 2
        #plt.plot(10*radii**(-1), p_av, 'o', mfc='none', label=folder + ', $\epsilon_r$ = {:1.2f}'.format(epsilon), color='black') # all points
        # style 3
        #plt.plot(10*radii**(-1), p_av, '.', label=folder + ', $\epsilon_r$ = {:1.2f}'.format(epsilon), color='black') # all points
        # style 4
        #plt.plot(10*radii**(-1), p_av, label=folder + ', $\epsilon_r$ = {:1.2f}'.format(epsilon), color='black') # all points
        #plt.plot(10*radii**(-1), p_av, label= '$\Delta E_{env} + \Delta V_{core} + \Delta V_{env}$', color='black') # all points
        plt.plot(10*radii**(-1), p_av, label= '$P_{env}$', color='black') # all points
        #plt.plot(10*radii_not_renormalized**(-1), p_av, label=folder + ', $\epsilon_r$ = {:1.2f}'.format(epsilon), color='grey') # not renormalized
        #plt.text(0.05, 1.6, folder, fontsize=12, bbox={'facecolor': 'white', 'edgecolor': 'white', 'alpha': 1, 'pad':10})


        eps_range = np.zeros(n)
        for i in range(n):
            eps_range[i] = eps_only(folder, radii, p_av, [lower_limits_eps[i], upper_limits_eps[i]])
        eps_std = np.std(eps_range)
        eps_mean = np.mean(eps_range)

        print('eps_mean ' + folder + '= ', eps_mean)
        print('eps_std ' + folder + '= ', eps_std)

        plt.xlim([0, 2])
        plt.ylim([-1, 2])
        plt.grid()
        plt.ylabel('Components of polarization energy, eV')
        plt.xlabel('$10/R,  10 / \AA^{-1} $')
        if radii_renormalization:
            plt.xlabel('Inverse distance renormalized $10/R,  10 / \AA^{-1} $')
        else:
            plt.xlabel('Inverse distance $10/R, \AA^{-1} $')
        plt.legend()
        ax1 = plt.gca()
        add_inverse_axis(ax1, rs_plot=np.array([5, 7, 10, 20, 30, 50]))
        plt.text(0.1, 1.6, folder, fontsize=12, bbox={'facecolor': 'white', 'edgecolor': 'grey', 'alpha': 0.5, 'pad':10})
        plt.savefig('p_average_componentwise_{}.png'.format(folder), dpi=600, bbox_inches='tight')



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