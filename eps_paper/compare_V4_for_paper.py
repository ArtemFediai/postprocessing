"""
NEWER VERSION!
Compute/plot:
IP(R)
EA(R)
0.5(IP(R) - EA(R))
as extracted from QP (damped shells)
"""

from __future__ import print_function
from __future__ import absolute_import

from typing import List

from IP_or_EA_vs_R_ver_4 import QPOutput, Timer
from IP_or_EA_vs_R_ver_4 import add_inverse_axis
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
from palettable.colorbrewer.sequential import Blues_5 as anion_palete
from palettable.colorbrewer.sequential import Reds_5 as cation_palete
from palettable.colorbrewer.sequential import Greys_5 as black_palete


def main():
    with Timer() as t:

        # for eps
        ab = [20,30]
        #

        folders = os.listdir('.')
        folders: List[str] = [folder for folder in folders if os.path.exists(folder+'/Analysis') and os.path.exists(folder+'/Analysis'+'/IP_by_radius')]
#        A,B,C = [],[],[]
        # compare_dict = {'dE0_cation': [], 'dE0_anion': [],
        #                 'dV0_cation': [], 'dV0_anion': [],
        #                 'dVe_cation': [], 'dVe_anion': [],
        #                 'dEe_cation': [], 'dEe_anion': [],
        #                 'list_radii': [], 'list_folders': []}
        quick_dir = {}

        plt.figure()
        for i, f in enumerate(folders):
            os.chdir(f)
            qp_settings_file = 'settings_ng.yml'  # default name of the QP settings file
            with open(qp_settings_file, 'r') as fid:
                qp_settings = yaml.load(fid, Loader=yaml.FullLoader)
            R = int(qp_settings['System']['Shells']['0']['cutoff'])
            #ab = [R // 2, R]  # ab is the left/right limit of the epsilon evaluation interval
#            my_mol_name = os.listdir(f + '/Analysis')[0].split('_')[1]  # will work only for 1-component material #TODO: wrong

            # here you can re-define R, ab, my_mol_name if necessary -->
            # R = 40
            # ab = [20, 40]
            # my_mol_name ='0ea4d81ac970d3f4fdbbe46acd91a041'
            # <--

            # get ea/ip output
            Analysis_output = QPOutput()  # ini: returns target folders

            Analysis_output.extract_energies()

            Analysis_output.extract_true_vacuum_energies()

            Analysis_output.get_polarization_elementwise()

            Analysis_output.plot_p()

            Analysis_output.plot_IPEA()

            Analysis_output.save_polarization_elementwise()

            Analysis_output.extract_eps_from_polarization(ab=ab)

            os.chdir('../')

            #plt.plot(1./Analysis_output.inv_radii*10, Analysis_output.dEe_anion, label = f)

            quick_dir[f] = Analysis_output

            # A.append(Analysis_output.dEe_anion)
            # B.append(Analysis_output.dVe_anion)
            # C.append(Analysis_output.dV0_anion)

#
        xlim = [0, 2]

        for i, f in enumerate(folders):
            r = quick_dir[f].inv_radii * 10
            x = quick_dir[f]
            sum = x.dEe_anion + x.dVe_anion + x.dV0_anion + x.dE0_anion
            plt.plot(r, x.dEe_anion, color='C'+str(i))
            plt.plot(r, x.dVe_anion, color='C'+str(i))
            plt.plot(r, x.dV0_anion, color='C'+str(i))
            plt.plot(r, x.dE0_anion*np.ones(len(r)), color='C'+str(i))
            plt.plot(r, sum, linewidth='1.5', color='C'+str(i), label = f)

        plt.legend(folders)

        plt.ylabel('dEe')
        plt.xlabel('1/r')
        plt.legend()
        plt.xlim(xlim)
        #plt.grid()
        ax1 = plt.gca()
        add_inverse_axis(initial_axis=ax1)
        #plt.grid()
        plt.savefig('P_componentswise.png', dpi=600)
        plt.close()


        # Figure 2
        plt.figure()
        for i, f in enumerate(folders):
            x = quick_dir[f]
            r = x.inv_radii*10
            plt.plot(r,0.5*(x.P_cation+x.P_anion), label = f + '; eps = {:.2f}'.format((x.epsilon)))

        plt.legend(folders)
        plt.ylabel('$0.5(P^{+} + P^{-})$')
        plt.xlabel('1/r, $\AA^{-1}$')
        plt.legend()
        plt.xlim(xlim)
        #plt.grid()
        ax1 = plt.gca()
        add_inverse_axis(initial_axis=ax1)
        #plt.grid()
        plt.savefig('P_average.png', dpi=600)
        plt.close()


        # Figure 3 P+ P- total guys
        plt.figure()

        ##
        ##

        for i, f in enumerate(folders):
            x = quick_dir[f]
            r = x.inv_radii*10
            plt.plot(r,x.P_cation, label = f + ' (cation)')
            plt.plot(r,x.P_cation, label = f + ' (anion)')

        plt.legend(folders)
        plt.ylabel('$0.5(P^{+} + P^{-})$')
        plt.xlabel('1/r, $\AA^{-1}$')
        plt.legend()
        plt.xlim(xlim)
        #plt.grid()
        ax1 = plt.gca()
        add_inverse_axis(initial_axis=ax1)
        #plt.grid()
        plt.savefig('P_for_paper.png', dpi=600)
        plt.close()


        plt.figure(figsize=[4.5, 3.5])
        xlim = [0, 3]
        for i, f in enumerate(reversed(folders)):
            x = quick_dir[f]
            r = x.inv_radii * 10
            #
            IP_vac = x.E0_plus_vacuum['mean'] - x.E0_0_vacuum['mean']
            IP_0 = x.E0_plus['mean'] - x.E0_0['mean']
            ddE0_plus = np.ones(len(x.P_cation)) * (IP_vac - IP_0)

            EA_vac = -x.E0_minus_vacuum['mean'] + x.E0_0_vacuum['mean']
            EA_0 = -x.E0_minus['mean'] + x.E0_0['mean']
            ddE0_minus = np.ones(len(x.P_anion)) * (EA_vac - EA_0)
            #
            if len(folders)-1 != i:
                plt.plot(r, x.P_anion - x.dE0_anion + ddE0_minus*0, color=anion_palete.mpl_colors[i+1])
                plt.plot(r, x.P_cation - x.dE0_cation + ddE0_plus*0, color=cation_palete.mpl_colors[i+1])
                plt.plot(r, 0.5*(x.P_cation - x.dE0_cation + x.P_anion - x.dE0_anion), color=black_palete.mpl_colors[i+1])
            else:
                plt.plot(r, x.P_anion -x.dE0_anion + ddE0_minus*0, color=anion_palete.mpl_colors[i+1], label='$\mathrm{P}^{-}$')
                plt.plot(r, x.P_cation - x.dE0_cation + ddE0_plus*0, color=cation_palete.mpl_colors[i+1], label='$\mathrm{P}^{+}$')
                plt.plot(r, 0.5*(x.P_cation - x.dE0_cation + x.P_anion -x.dE0_anion), color=black_palete.mpl_colors[i+1], label='$\mathrm{P}$')

        plt.ylabel('$\mathrm{P}^{+}, \mathrm{P}^{-}, \mathrm{P}$')
        plt.xlabel('$1/R,  10 / \AA^{-1} $')
        plt.legend()
        plt.xlim(xlim)
        ax1 = plt.gca()
        add_inverse_axis(initial_axis=ax1, rs_plot=np.array([10,20,30,40,50,60]), rs_grid=np.array([10,20,30,40,50,60]))
        plt.savefig('P_for_paper.png', dpi=600, bbox_inches='tight')
        plt.savefig('P_for_paper.svg', bbox_inches='tight')
        plt.close()


    print("\nI am done")
    print("Total Computation Time: {} sec".format(t.interval))


if __name__ == '__main__':
    main()
