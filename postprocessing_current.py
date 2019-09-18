"""
this is in fact traps simulations n_r n_t
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
from shutil import copyfile


def main():
    """
    "1"             READ
    """
    error = 'std'  # std_err or std

    n_r, n_t = [int(np.loadtxt('n_r_n_t.dat')[0]), int(np.loadtxt('n_r_n_t.dat')[1])]
    field = 0.036 * 10 ** 9  # field, V/nm SI
    n_traps = np.loadtxt('N_t.dat')
    #n_t = 19 #  hc
    #n_traps = n_traps[0:n_t]

    """
    "2"         COMPUTE AND SAVE: AVERAGE DISTRIBUTIONS
    """

    my_current = Current(n_t, n_r)  # ini with empties

    if not os.path.isdir('postprocessing'):
        os.makedirs('postprocessing')
    if not os.path.isfile("postprocessing/n_traps.dat"):
        np.savetxt("postprocessing/n_traps.dat", n_traps)  # ugly

    if (os.path.isfile('postprocessing/current.dat') and
        os.path.isfile('postprocessing/mobility.dat') and
        os.path.isfile('postprocessing/conductivity.dat')and
        os.path.isfile('postprocessing/std_err_cu.dat') and
        os.path.isfile('postprocessing/std_err_mo.dat') and
        os.path.isfile('postprocessing/std_err_co.dat')):

        my_current.current = np.loadtxt('postprocessing/current.dat')
        my_current.mobility = np.loadtxt('postprocessing/mobility.dat')
        my_current.conductivity = np.loadtxt('postprocessing/conductivity.dat')
        my_current.std_err_cu = np.loadtxt('postprocessing/std_err_cu.dat')
        my_current.std_err_mo = np.loadtxt('postprocessing/std_err_mo.dat')
        my_current.std_err_co = np.loadtxt('postprocessing/std_err_co.dat')
    else:
        my_current.compute(field)

    """
     3          PLOT
    """
    print("current: \n", my_current.current)

    plt.figure()
    plt.plot(n_traps, my_current.current)
    plt.errorbar(n_traps, my_current.current, yerr=my_current.std_err_cu, fmt='o-')
#    plt.xscale('log')
    plt.yscale('log')
    plt.ylabel('current, A/m$^2$')
    plt.xlabel('trap molar fraction')
    plt.savefig("current.png")
    plt.close()

    plt.figure()
    if error == 'std':
        plt.errorbar(n_traps, my_current.mobility*1E4, yerr=(my_current.std_err_mo*1E4*n_r)) #  std
    else:
        plt.errorbar(n_traps, my_current.mobility*1E4, yerr=my_current.std_err_mo*1E4, fmt='o-') # std_err
    plt.xscale('log')
    plt.yscale('log')
    plt.ylabel('mobility, cm$^2$V$^-1$s$^-1$')
    plt.xlabel('trap molar fraction')
    plt.savefig("mobility.png")
    plt.close()

    plt.figure()
    plt.errorbar(n_traps, my_current.conductivity, yerr=my_current.std_err_co, fmt='o-')
    plt.ylabel('Conductivity, S/cm')
    plt.xscale('log')
    plt.yscale('log')
    plt.ylabel('conductivity, S/m')
    plt.xlabel('trap molar fraction')
    plt.savefig("conductivity.png")
    plt.close()

    """
    N-1     PRINT
    """

    # uncomment to print

    print("current: \n {} \n".format(my_current.current))
    print("mobility: \n {} \n".format(my_current.mobility))
    print("conductivity: \n {} \n".format(my_current.conductivity))

    """
    N       PLOT
    """


class Current:
    def __init__(self, n_t, n_r):
        self.n_t = n_t
        self.n_r = n_r
        self.current = np.zeros(n_t)
        self.mobility = np.zeros(n_t)
        self.conductivity = np.zeros(n_t)
        self.std_err_cu = np.zeros(n_t)
        self.std_err_mo = np.zeros(n_t)
        self.std_err_co = np.zeros(n_t)

    def compute(self, f):
        fn = "t_{}/r_{}/output_job_0"
        current, mobility, field, conductivity = np.zeros(self.n_r), np.zeros(self.n_r), np.zeros(self.n_r), np.zeros(self.n_r)
        for i_t in range(0, self.n_t):
            for i_r in range(0, self.n_r):
                print("read-in i_t = {}, i_r = {}".format(i_t, i_r))
                current[i_r] = extract.extract_e(fn.format(i_t, i_r), "final current density: ")
                mobility[i_r] = extract.extract_e(fn.format(i_t, i_r), "final mobility: ")
                field[i_r] = f
                conductivity[i_r] = current[i_r]/field[i_r]
            self.current[i_t] = gmean(current)
            self.mobility[i_t] = gmean(mobility)
            self.conductivity[i_t] = gmean(conductivity)

            #self.current[i_t] = np.mean(current)
            #self.mobility[i_t] = np.mean(mobility)g
            #self.conductivity[i_t] = np.mean(conductivity)
            self.std_err_cu[i_t] = np.std(current)/self.n_r
            self.std_err_mo[i_t] = np.std(mobility)/self.n_r
            self.std_err_co[i_t] = np.std(conductivity)/self.n_r

        np.savetxt('postprocessing/current.dat', self.current)
        np.savetxt('postprocessing/mobility.dat', self.mobility)
        np.savetxt('postprocessing/conductivity.dat', self.conductivity)
        np.savetxt('postprocessing/std_err_cu.dat', self.std_err_cu)
        np.savetxt('postprocessing/std_err_mo.dat', self.std_err_mo)
        np.savetxt('postprocessing/std_err_co.dat', self.std_err_co)


if __name__ == '__main__':
    main()
