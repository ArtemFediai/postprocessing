# ad hoc
# needs the existing of mentioned files

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

pcs_pcs_folder = 'postprocessing_output/pcs_pcs'
pcs_qp_folder = 'postprocessing_output/pcs_qp'

#for aNPD

materials = ['aNPD', 'C60', 'TCTA']

for i,material in enumerate(materials):
    plt.figure(material)
    p_av_pcs_pcs = np.loadtxt(pcs_pcs_folder + '/' + 'p_av_' + material + '.dat')
    radii_pcs_pcs = np.loadtxt(pcs_pcs_folder + '/' + 'radii.dat')
    p_av_pcs_qp = np.loadtxt(pcs_qp_folder + '/' + 'p_av_' + material + '.dat')
    radii_pcs_qp = np.loadtxt(pcs_qp_folder + '/' + 'radii.dat')

    plt.plot(10/radii_pcs_pcs, p_av_pcs_pcs, label='pcs-pcs', marker='o')
    plt.plot(10/radii_pcs_qp, p_av_pcs_qp, label='pcs-DFT', marker='x')

    plt.xlim([0, 2])
    plt.ylim([0, 1.0])
    plt.ylim(bottom=0)
    plt.grid()
    plt.ylabel('Average polarization energy, eV')
    plt.xlabel('$10/R,  10 / \AA^{-1} $')
    # if radii_renormalization:
    #     plt.xlabel('Inverse distance renormalized $1/R,  10 / \AA^{-1} $')
    # else:
    #     plt.xlabel('Inverse distance $1/R,  10 / \AA^{-1} $')
    plt.legend()
    ax1 = plt.gca()
    add_inverse_axis(ax1)

    plt.savefig('postprocessing_output/' + 'p_av_compare_{}.png'.format(material), dpi=600, bbox_inches='tight')

    print('different in eV {} for {}'.format((p_av_pcs_qp[-1] - p_av_pcs_pcs[-1]), material))
