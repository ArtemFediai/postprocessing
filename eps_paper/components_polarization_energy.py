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
import seaborn as sns; sns.set()
import pandas as pd
#sns.set_theme(style="whitegrid")
sns.set_style('whitegrid')
sns.color_palette()

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

    mean_array_of_ps_at_0 = np.array([0.0, 0.0, 0.0])
    store_it = np.zeros([3,3])
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
        delta_delta = 0.5 * (delta_delta_cation - delta_delta_anion)

        p_av = (e0i_anion + eii_anion + v0i_anion + e0i_cation + eii_cation + v0i_cation)/2.0
        deii = 0.5*(e0i_anion+e0i_cation) #dE_env
        dv0i = 0.5*(v0i_anion+v0i_cation) #dR_core
        dvii = 0.5*(eii_anion+eii_cation) #dV_env
        de0 = 0.5*(e0_anion + e0_cation)

        bulks = []

        #total ->
        point_a, point_b, point_b1, target_i, point_a1, epsilon = eps_mod(folder, radii, radii_not_renormalized, p_av, [20, 50], [1, 1E5])
        plt.plot([point_a[0], point_b[0]], [point_a[1], point_b[1]], LineStyle='dotted', color='black')
        plt.plot(point_a[0], point_a[1], LineStyle='dotted', marker='o', color='black')
        #<- total

        #partial ->
        array_of_ps_at_0 = []
        array_of_ps = [deii, dv0i, dvii] # various polarizations
        for j, p in enumerate(array_of_ps):
            point_a, point_b, point_b1, _, point_a1, _ = eps_mod(folder, radii, radii_not_renormalized, p, [20, 40], [1, 1E5])
            print([point_a1[1]])
            array_of_ps_at_0.append(point_a1[1])
            plt.plot([point_a[0], point_b[0]], [point_a[1], point_b[1]], LineStyle='dotted', color=f'C{j}')
            plt.plot([point_a[0]], [point_a[1]], LineStyle='dashed', marker='o', color=f'C{j}')
            bulks.append(point_a[1])

        bulks.append(delta_delta)##dd

        print(array_of_ps_at_0)
        #<- partial

        plt.plot(10*radii**(-1), deii, label='$\Delta E_\mathrm{env}$', color='C0')#
        plt.plot(10*radii**(-1), dv0i, label='$\Delta V_\mathrm{core}$', color='C1') #
        plt.plot(10*radii**(-1), dvii, label='$\Delta V_\mathrm{env}$', color='C2') #
        plt.plot(10*radii**(-1), p_av, label= '$\mathrm{P}_\mathrm{env}$', color='black')  # all points

        eps_range = np.zeros(n)
        for j in range(n):
            eps_range[j] = eps_only(folder, radii, p_av, [lower_limits_eps[j], upper_limits_eps[j]])
        eps_std = np.std(eps_range)
        eps_mean = np.mean(eps_range)

        print('eps_mean ' + folder + '= ', eps_mean)
        print('eps_std ' + folder + '= ', eps_std)

        plt.xlim([-0.1, 2])
        plt.ylim([-1, 2.5])
        plt.ylabel('Components of $\mathrm{P}$, eV')
        plt.xlabel('$10/R,  10 / mathrm{\AA}^{-1} $')
        if radii_renormalization:
            plt.xlabel('Inverse renormalized radius $10/R,  10 / \mathrm{\AA}^{-1} $')
        else:
            plt.xlabel('Inverse polarization radius $10/R, \mathrm{\AA}^{-1} $')
        plt.legend()
        ax1 = plt.gca()
        add_inverse_axis(ax1, rs_plot=np.array([5, 7, 10, 20, 30, 50]), x_label='$R, \mathrm{\AA}$')

        plt.grid()
        plt.text(0.6, 2.15, folder, fontsize=12, bbox={'facecolor': 'None', 'edgecolor': 'None', 'alpha': 0.5, 'pad':10})
        plt.savefig('p_average_componentwise_{}.png'.format(folder), dpi=600, bbox_inches='tight')
        plt.savefig('p_average_componentwise_{}.svg'.format(folder), bbox_inches='tight')
        plt.close()


        #<-- pie chart
        data = {'$\Delta E_\mathrm{env}$': abs(bulks[0]),
                '$\Delta V_\mathrm{core}$': abs(bulks[1]),
                '$\Delta V_\mathrm{env}$': abs(bulks[2]),
                '$\Delta\Delta E_0$': abs(bulks[3])}

        # data = {'$\Delta E_\mathrm{env}$': abs(bulks[0]),
        #         '$\Delta V_\mathrm{core}$': abs(bulks[1]),
        #         '$\Delta V_\mathrm{env}$': abs(bulks[2])}

        pie, ax = plt.subplots(figsize=[3.5, 4.5])
        labels = data.keys()
        plt.pie(x=list(data.values()), autopct="%.1f%%", explode=[0.05] * len(data), labels=labels, pctdistance=0.5)
        plt.title(folder, fontsize=14)
        pie.savefig(f"PieChart_{folder}.png", dpi=600)
        pie.savefig(f"PieChart_{folder}.svg")
        plt.close()

        #--> bar chart
        data_signed = \
                {'$\Delta V_\mathrm{core}$': (bulks[1]),
                '$+\Delta E_\mathrm{env}$': (bulks[0]),
                '$+\Delta V_\mathrm{env}$': (bulks[2]),
                '$+\Delta\Delta E_0$': (bulks[3])}

        data_cumulative = {}
        data_numpy = np.array([x for x in data_signed.values()])
        data_str = [x for x in data_signed.keys()]
        data_cumulative = data_numpy.cumsum()

        # data_for_panda = [[x,y] for x,y in zip(data_str, data_numpy)]
        data_for_panda = [[x,y] for x,y in zip(data_str, data_cumulative)]
        df = pd.DataFrame(data=data_for_panda, columns=["Cumulative sum", "Polarization"])
        ax = sns.barplot(data=df, x='Cumulative sum', y='Polarization', palette='Blues_d')
        ax.set(ylabel="Cumulative polarization, eV")
        ax.set(xlabel="Polarization components")
        ax.set(ylim=[0, 3])
        plt.title(folder, fontsize=14)
        plt.savefig(f'bar_plot_{folder}.png', dpi=600)
        plt.savefig(f'bar_plot_{folder}.svg')
        plt.close()
        #<-- bar chart
        normalized = np.abs(np.array(array_of_ps_at_0))  / np.sum(np.abs(array_of_ps_at_0))
        mean_array_of_ps_at_0 +=  normalized

        store_it[i,:] = normalized


        # std for components

    print("computing mean")
    print(mean_array_of_ps_at_0/len(folders))
    print("all")
    print(store_it)
    for ii in range(np.shape(store_it)[1]):
        print(store_it[:, ii])
        print(np.std(store_it[:,ii]))


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