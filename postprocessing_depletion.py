"""
this is for simulations n_dop n_t. Extracting depletion parameters and mobility
"""

from __future__ import print_function
from __future__ import absolute_import
import extract
import numpy as np
import matplotlib

matplotlib.use('pdf')
import matplotlib.pyplot as plt
import os
from scipy.stats.mstats import gmean as gmean
from extract import num_fs as num_fs


def main():
    """
    "0"           SHORT MANUAL INPUT
    """
    folders_patterns = ['dop_', 'r_']  # array of strings i.e. ['dop_', 'r_']
    error = 'std_err'  # std_err or std

    """
    "1"             READ
    """
    n_d, n_r = num_fs(folders_patterns)
    print("n_d, n_r = {}, {}".format(n_d, n_r))

    dop = np.loadtxt('doping.txt')  # all doping

    #n_d = 6
    #dop = dop[0:n_d]
    #n_r = 10
    """
    "2"         COMPUTE AND SAVE: AVERAGE DISTRIBUTIONS
    """

    my_depletion = Depletion(n_d, n_r)  # ini with empties

    if not os.path.isdir('postprocessing'):
        os.makedirs('postprocessing')
    if not os.path.isfile("postprocessing/dop.dat"):
        np.savetxt("postprocessing/dop.dat", dop)  # ugly

    if (os.path.isfile('postprocessing/depl_length.dat') and
            os.path.isfile('postprocessing/mobile_fraction.dat') and
            os.path.isfile('postprocessing/std_err_de.dat') and
            os.path.isfile('postprocessing/std_err_mo.dat')):

        my_depletion.depl_length = np.loadtxt('postprocessing/depl_length.dat')
        my_depletion.mobile_fraction = np.loadtxt('postprocessing/mobile_fraction.dat')
        my_depletion.std_err_de = np.loadtxt('postprocessing/std_err_de.dat')
        my_depletion.std_err_mo = np.loadtxt('postprocessing/std_err_mo.dat')
    else:
        my_depletion.compute()

    """
     3          PLOT
    """
    # Fig 1: dop vs. L | normal scales
    plt.figure()
    plt.plot(dop, my_depletion.depl_length)
    plt.errorbar(dop, my_depletion.depl_length, yerr=my_depletion.std_err_de, fmt='o-')
    #plt.xscale('log')
    plt.ylabel('Depletion length, nm')
    plt.xlabel('Doping')
    plt.savefig("Depletion_length.png")
    plt.close()

    # include vacuum
    w = my_depletion.depl_length
    delta = 0.8
    coeff = w**2/(delta+w)**2
    #

    # Fig 2: dop. vs. mobile fraction | xscale log
    plt.figure()
    if error == 'std':
        plt.errorbar(dop, my_depletion.mobile_fraction, yerr=(my_depletion.std_err_mo * n_r), capsize=3)  # std
    else:
        plt.errorbar(dop, my_depletion.mobile_fraction, yerr=my_depletion.std_err_mo, capsize=3, label='as extracted')  # std_err
        plt.errorbar(dop, my_depletion.mobile_fraction*coeff, yerr=my_depletion.std_err_mo, capsize=3, label='corrected')  # std_err
    plt.xscale('log')
    plt.ylim(bottom=0)
    plt.ylabel('Fraction of mobile charge carriers')
    plt.xlabel('Doping')
    plt.legend()
    plt.savefig("MobileFraction.png")
    plt.close()

    # Fig. 3: dop**(-1/2) vs L
    plt.figure()
    #plt.plot(dop, my_depletion.depl_length)
    plt.errorbar((dop)**(-1/2), my_depletion.depl_length, yerr=my_depletion.std_err_de, fmt='o')
    plot_trend((dop)**(-1/2), my_depletion.depl_length, color_num='0')
    plt.ylabel('Depletion length, nm')
    plt.xlabel('Dopant molar fraction$^{-1/2}$')
    plt.savefig("Depletion_length_vs_sqrt_doping.png")
    plt.close()



    """
    N-1     PRINT
    """

    # uncomment to print

    print("Depletion: \n {} \n".format(my_depletion.depl_length))
    print("Mobile fraction: \n {} \n".format(my_depletion.mobile_fraction))

    """
    N       PLOT
    """


class Depletion:
    def __init__(self, n_d, n_r):
        self.n_d = n_d
        self.n_r = n_r
        self.depl_length = np.zeros(n_d)
        self.mobile_fraction = np.zeros(n_d)
        self.std_err_de = np.zeros(n_d)
        self.std_err_mo = np.zeros(n_d)

    def compute(self):
        fn = "dop_{}/r_{}/output_job_0"  # ugly
        depl_length, mobile_fraction = np.zeros(self.n_r), np.zeros(self.n_r)
        for i_d in range(0, self.n_d):
            for i_r in range(0, self.n_r):
                print("read-in i_dop = {}, i_r = {}".format(i_d, i_r))
                depl_length[i_r] = extract.extract_float(fn.format(i_d, i_r), "depletion length: ")
                mobile_fraction[i_r] = extract.extract_float(fn.format(i_d, i_r), "mobile carrier fraction: ")
            #self.depl_length[i_d] = gmean(-depl_length)
            #self.mobile_fraction[i_d] = gmean(mobile_fraction)

            self.depl_length[i_d] = np.mean(-depl_length)
            self.mobile_fraction[i_d] = np.mean(mobile_fraction)
            self.std_err_de[i_d] = np.std(depl_length) / self.n_r
            self.std_err_mo[i_d] = np.std(mobile_fraction) / self.n_r

        np.savetxt('postprocessing/depl_length.dat', self.depl_length)
        np.savetxt('postprocessing/mobile_fraction.dat', self.mobile_fraction)
        np.savetxt('postprocessing/std_err_de.dat', self.std_err_de)
        np.savetxt('postprocessing/std_err_mo.dat', self.std_err_mo)

def plot_trend(x,y,color_num = '0'):
    z = np.polyfit(x, y, 1)
    p = np.poly1d(z)
    plt.plot(x, p(x), color='C{}'.format(color_num), linestyle=':')


if __name__ == '__main__':
    main()
