"""
input format:
i_t/i_r that is the trap number and the replica number
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
    TrapLevelExists = False
    A = np.loadtxt('n_r_n_dop.txt')  # num of replicas, num of dopings
    n_r, n_d = [int(np.loadtxt('n_r_n_dop.txt')[0]), int(np.loadtxt('n_r_n_dop.txt')[1])]
    # hard coded
    scale = 0.196  # this is only once for the paper
    # end: hard-coded
    field = 0.04 * 10 ** 9  # field, V/nm SI
    dop = np.loadtxt('doping.txt')


    """
    "2"         COMPUTE AND SAVE: AVERAGE DISTRIBUTIONS
    """

    my_current = Current(n_d, n_r)  # ini with empties

    if not os.path.isdir('postprocessing'):
        os.makedirs('postprocessing')
    if not os.path.isfile("postprocessing/doping.dat"):
        np.savetxt("postprocessing/doping.dat", dop)  # ungly

    if (os.path.isfile('postprocessing/current.dat') and
        os.path.isfile('postprocessing/mobility.dat') and
        os.path.isfile('postprocessing/conductivity.dat')and
        os.path.isfile('postprocessing/std_err_cu.dat') and
        os.path.isfile('postprocessing/std_err_mo.dat') and
        os.path.isfile('postprocessing/std_err_co.dat')):

        my_current.current = np.loadtxt('postprocessing/current.dat')
        my_current.mobility = np.loadtxt('postprocessing/mobility.dat')
        my_current.conductivity = np.loadtxt('postprocessing/conductivity.dat')
        my_current.std_err_current = np.loadtxt('postprocessing/std_err_cu.dat')
        my_current.std_err_mobility = np.loadtxt('postprocessing/std_err_mo.dat')
        my_current.std_err_conductivity = np.loadtxt('postprocessing/std_err_co.dat')
    else:
        my_current.compute(field)

    """
     3          PLOT
    """
    print("current: \n", my_current.current)

    plt.figure()
    plt.plot(dop, my_current.current*scale)
    plt.xscale('log')
    plt.yscale('log')
    plt.savefig("current.png")

    plt.figure()
    plt.plot(dop, my_current.mobility*scale)
    plt.xscale('log')
    plt.yscale('log')
    plt.savefig("mobility.png")

    plt.figure()
    plt.plot(dop, my_current.conductivity*scale*100)
    plt.ylabel('Conductivity, S/cm')
    plt.xscale('log')
    plt.yscale('log')
    plt.savefig("conductivity.png")

    """
    N-1     PRINT
    """

    # uncomment to print
    """
    print("LUMO mean", level_lumo_plus.mean)
    """
    """
    print("current: \n {} \n".format(MyCurrent.current))
    print("mobility: \n {} \n".format(MyCurrent.mobility))
    print("conductivity: \n {} \n".format(MyCurrent.conductivity))
    """
    """
    N       PLOT
    """


class Current:
    def __init__(self, n_d, n_r):
        self.n_d = n_d
        self.n_r = n_r
        self.current = np.zeros(n_d)
        self.mobility = np.zeros(n_d)
        self.conductivity = np.zeros(n_d)
        self.std_err_cu = np.zeros(n_d)
        self.std_err_mo = np.zeros(n_d)
        self.std_err_co = np.zeros(n_d)

    def compute(self, f):
        fn = "dop_{}/r_{}/output_job_0"
        current, mobility, field, conductivity = np.zeros(self.n_r), np.zeros(self.n_r), np.zeros(self.n_r), np.zeros(self.n_r)
        for i_d in range(0, self.n_d):
            for i_r in range(0, self.n_r):
                current[i_r] = extract.extract_e(fn.format(i_d, i_r), "final current density: ")
                mobility[i_r] = extract.extract_e(fn.format(i_d, i_r), "final mobility: ")
                #field[i_r] = extract.extract_float(fn.format(i_d, i_r),"applying field of [ ")  #in V/nm
                field[i_r] = f
                conductivity[i_r] = current[i_r]/field[i_r]
            self.current[i_d] = gmean(current)
            self.mobility[i_d] = gmean(mobility)
            self.conductivity[i_d] = gmean(conductivity)

            #self.current[i_d] = np.mean(current)
            #self.mobility[i_d] = np.mean(mobility)
            #self.conductivity[i_d] = np.mean(conductivity)
            self.std_err_cu[i_d] = np.std(current)/self.n_r
            self.std_err_mo[i_d] = np.std(mobility)/self.n_r
            self.std_err_co[i_d] = np.std(conductivity)/self.n_r

        np.savetxt('postprocessing/current.dat', self.current)
        np.savetxt('postprocessing/mobility.dat', self.mobility)
        np.savetxt('postprocessing/conductivity.dat', self.conductivity)
        np.savetxt('postprocessing/std_err_cu.dat', self.std_err_cu)
        np.savetxt('postprocessing/std_err_mo.dat', self.std_err_mo)
        np.savetxt('postprocessing/std_err_co.dat', self.std_err_co)


if __name__ == '__main__':
    main()
