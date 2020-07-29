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
from IP_or_EA_vs_R_ver_4 import QPOutput, Timer
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


def main():
    with Timer() as t:
        folders = os.listdir('.')
        folders = [folder for folder in folders if os.path.exists(folder+'/Analysis')]
        A,B,C = [],[],[]
        plt.figure()
        for i, f in enumerate(folders):
            os.chdir(f)
            qp_settings_file = 'settings_ng.yml'  # default name of the QP settings file
            with open(qp_settings_file, 'r') as fid:
                qp_settings = yaml.load(fid, Loader=yaml.FullLoader)
            R = int(qp_settings['System']['Shells']['0']['cutoff'])
            ab = [R // 2, R]  # ab is the left/right limit of the epsilon evaluation interval
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

            os.chdir('../')

            #plt.plot(1./Analysis_output.inv_radii*10, Analysis_output.dEe_anion, label = f)
            A.append(Analysis_output.dEe_anion)
            B.append(Analysis_output.dVe_anion)
            C.append(Analysis_output.dV0_anion)

#
        xlim = [0, 1.25]
        for i,f in enumerate(folders):
            plt.plot(Analysis_output.inv_radii*10, A[i], label=f)
        plt.ylabel('dEe')
        plt.xlabel('1/r')
        plt.legend()
        plt.xlim(xlim)
        plt.savefig('dEe.png', dpi=600)
        plt.close()
#
        for i,f in enumerate(folders):
            plt.plot(Analysis_output.inv_radii*10, B[i], label=f)
        plt.ylabel('dVe')
        plt.xlabel('1/r')
        plt.legend()
        plt.xlim(xlim)
        plt.savefig('dVe.png', dpi=600)
        plt.close()
#
        for i,f in enumerate(folders):
            plt.plot(Analysis_output.inv_radii*10, C[i], label=f)
        plt.ylabel('dV0')
        plt.xlabel('1/r')
        plt.legend()
        plt.xlim(xlim)
        plt.savefig('dV0.png', dpi=600)
        plt.close()

    print("\nI am done")
    print("Total Computation Time: {} sec".format(t.interval))


if __name__ == '__main__':
    main()
