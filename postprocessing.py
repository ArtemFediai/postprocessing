"""
1. Read
2. Compute average distributions
3. Get parameters of average distributions
4. Save distributions
"""
from __future__ import print_function
from __future__ import absolute_import
import extract
import numpy as np
import matplotlib
# matplotlib.use('pdf')
import matplotlib.pyplot as plt
import os
from shutil import copyfile


def main():
    """
    "1"
    """
    A = np.loadtxt('n_r_n_dop.txt')  # num of replicas, num of dopings
    n_r, n_d = [int(np.loadtxt('n_r_n_dop.txt')[0]), int(np.loadtxt('n_r_n_dop.txt')[1])]
    dop = np.loadtxt('doping.txt')
    a = np.loadtxt("a.txt")  # size of the system
    aaa = a*a*a

    N_r = np.linspace(0, n_r - 1, n_r, dtype=int)  # if custom N_r is required
    n_r = len(N_r)
    E1 = -7.0
    E2 = -3.0
    NE = 1000
    Energy = np.linspace(E1, E2, NE)  # this must be saved later!
    T = 300

    # ini
    HOMO = np.zeros([n_d, NE], dtype=np.float)
    LUMO_plus = np.zeros([n_d, NE], dtype=np.float)
    # end ini

    """
    "2"
    """
    if not os.path.isdir('postprocessing'):
        os.makedirs('postprocessing')
        
    if os.path.isfile('postprocessing/HOMO.dat'):
        HOMO[:, :] = np.loadtxt('postprocessing/HOMO.dat')
    else:
        HOMO[:, :] = get_average_distribution('IP_0', n_d, n_r, NE)  # this is an average distribution of HOMO
        for i_d in range(0, n_d):
            HOMO[i_d, :] /= aaa[i_d]
        np.savetxt('postprocessing/HOMO.dat', HOMO)

    if os.path.isfile('postprocessing/LUMO.dat'):
        LUMO_plus[:, :] = np.loadtxt('postprocessing/LUMO.dat')
    else:
        LUMO_plus[:, :] = get_average_distribution('EA_plus_0', n_d, n_r, NE)  # the same for LUMO+
        for i_d in range(0, n_d):
            LUMO_plus[i_d, :] /= aaa[i_d]
        np.savetxt('postprocessing/LUMO.dat', LUMO_plus)

    #quick check

    #end: quick check

    """
    "3"
    """

    level_homo = Level(HOMO, Energy)
    level_homo.compute_parameters()

    level_lumo_plus = Level(LUMO_plus, Energy)
    level_lumo_plus.compute_parameters()

    level_homo_plus_lumo = Level(HOMO + LUMO_plus, Energy)
    level_homo_plus_lumo.compute_parameters()

    """
    "4"
    """

    np.savetxt("postprocessing/HOMO_disorder.dat", level_homo.disorder)
    np.savetxt("postprocessing/LUMO_plus_disorder.dat", level_lumo_plus.disorder)
    np.savetxt("postprocessing/HOMO_plus_LUMO_disorder.dat", level_homo_plus_lumo.disorder)
    np.savetxt("postprocessing/energy.dat", Energy)
    np.savetxt("postprocessing/doping.dat", dop)


    """
    "5"         FERMI LEVEL
    """

    density_of_states = DOS(HOMO + LUMO_plus, Energy, T, dop)

    density_of_states.compute_parameters()
    density_of_states.compute_ef()

    print("ef = ", density_of_states.ef)
    print("LUMO mean", level_lumo_plus.mean)
    dE = Energy[-1]-Energy[-2]
    print("must be one:", np.sum(level_homo_plus_lumo.distribution[12, :])+dop[12])

    """
    OPTIONAL PLOTS
    """

    # uncomment to plot HOMO, LUMO and HOMO+LUMO disorder
    """
    plt.figure()
    plt.plot(dop, level_homo.disorder, label="HOMO")
    plt.plot(dop, level_lumo_plus.disorder, label="LUMO")
    plt.plot(dop, level_homo_plus_lumo.disorder, label="HOMO+LUMO")
    plt.xscale("log")
    plt.legend()
    # plt.show()
    plt.savefig("disorder.png")
    plt.close()
    # end: uncomment
    """

    # uncomment to plot HOMO and LUMO_plus
    """
    plt.figure()
    for i_d in range(n_d):
        plt.plot(Energy, HOMO[i_d, :])
    for i_d in range(n_d):
        plt.plot(Energy, LUMO_plus[i_d, :])
    plt.yscale('log')
    plt.show()
    exit()
    """
    # end: uncomment


class Level:
    def __init__(self, distribution, energy):
        self.distribution = distribution
        self.energy = energy
        self.n_d = np.size(distribution, 0)
        self.n_e = np.size(distribution, 1)
        self.disorder = np.zeros(self.n_d, dtype=np.float)
        self.max = np.zeros(self.n_d, dtype=np.float)
        self.mean = np.zeros(self.n_d, dtype=np.float)
        self.ef = np.zeros(self.n_d, dtype=np.float)

    def compute_parameters(self):
        for i_d in range(0, self.n_d):
            a = self.distribution[i_d, :]
            # self.disorder[i_d] = np.sqrt(np.cov(self.energy, aweights=a))
            self.max[i_d] = self.energy[a.argmax()]
            self.mean[i_d], self.disorder[i_d] = weighted_avg_and_std(self.energy, a)


class DOS(Level):
    """
    DOS must be sum of all relevant energy levels (HOMO and LUMO+)
    """
    def __init__(self, distribution, energy, temperature, doping):
        Level.__init__(self, distribution, energy)
        self.temperature = temperature
        self.doping = doping
        self.ef = np.zeros(self.n_d)

    def compute_ef(self):
        for i_d in range(0, self.n_d):
            self.ef[i_d] = self.ef_scalar(i_d)  # for a given doping compute Fermi level

    def ef_scalar(self, i_d):
        min_charge = 1E100  # big number
        i_e_min = 0  # not so nice
        for i_e in range(0, self.n_e):  # screen pot ef
            i1_ = self.I1(i_d)
            i2_ = self.I2(i_d, i_e)
            tmp = np.abs(i1_-i2_)
            #print(tmp)
            if tmp < min_charge:
                min_charge = tmp
                i_e_min = i_e  # index of the energy that minimize charge
        return self.energy[i_e_min]

    def I1(self, i_d):
        return self.doping[i_d]

    def I2(self, i_d, i_e):
        #de = self.energy[-1] - self.energy[-2]
        de = 1  # normalized in a way that the sum over all energies gives fraction of hosts ~< 1.0
        return np.sum(self.distribution[i_d, :] * (1.0-Fermi(self.energy, self.energy[i_e], self.temperature))) * de

def get_average_distribution(level_name, n_d, n_r, ne):
    # function to compute the average distribution for each doping
    # IP,EA--> [],_plus,_minus -->_X_0.dat   i.e. EA_plus_2_0.dat. Input: EA_plus
    tmp = np.zeros([n_d, ne])
    for i_d in range(0, n_d):
        for i_r in range(0, n_r):
            file_name = 'dop_{}/r_{}/experiments/tmp/{}_0.dat'.format(i_d, i_r, level_name)
            tmp[i_d, :] += np.loadtxt(file_name)
    return tmp / n_r


def weighted_avg_and_std(values, weights):
    average = np.average(values, weights=weights)
    variance = np.average((values - average) ** 2, weights=weights)
    return average, np.sqrt(variance)

def Fermi(e, EF, T):
    kB = 8.61733E-5#eV/K
    e = np.array(e)
    a = (e-EF)/(kB*T)
    return 1.0/(1.0+np.exp(a))

if __name__ == '__main__':
    main()
