"""
this is for simulations n_dop n_t. Extracting depletion parameters and mobility
"""

from __future__ import print_function
from __future__ import absolute_import
import extract
import numpy as np
import matplotlib
from scipy import constants as c
from ruamel.yaml import YAML

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
    dop = np.loadtxt('doping.txt')  # all doping. ugly. must be make if no postprocesing exists
    #n_d = 3
    #dop = dop[0:n_d-1]
    n_d = int(n_d)
    n_r = int(n_r)



    """
    "2"         COMPUTE AND SAVE: AVERAGE DISTRIBUTIONS
    """

    my_depletion = Depletion(n_d, n_r)  # ini Depletion obj with empties

    if not os.path.isdir('postprocessing'):
        os.makedirs('postprocessing')
    if os.path.isfile('postprocessing/depl.yaml'):
        fid = open('postprocessing/depl.yaml')
        yaml_obj = YAML(typ='safe')
        #from ruamel.yaml.comments import CommentedMap as OrderedDict
        #yaml_obj = YAML.load(fid)
        yaml_obj.register_class(Depletion)
        my_depletion = yaml_obj.load(fid)
        fid.close()
    else:
        my_depletion.compute()
        my_depletion.extract_yaml()
        yaml_obj = YAML(typ='safe')
        yaml_obj.register_class(Depletion)
        f2 = open('postprocessing/depl.yaml','w+')
        yaml_obj.dump(my_depletion, f2)
        f2.close()

    my_depletion.extract_level(level='IP_x', mol_num='0')


    """
    PLOT FROM YAML
    """

    # Fig 1: dop vs. L | normal scales





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
    coeff = np.array(w)**2/(np.array(delta)+np.array(w))**2
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

    # Fig. 4: thickness vs. L/2
    plt.figure()
    plt.plot(dop, np.array(my_depletion.thickness)/2, label='Half thickness')
    plt.plot(dop, my_depletion.depl_length, label='Depletion length')
    plt.ylabel('Depletion length, nm')
    plt.xlabel('Dopant molar fraction')
    plt.xscale('log')
    plt.legend()
    plt.savefig("Depletion_length_vs_Thickness.png")
    plt.close()

    # Fig. 5: thickness vs. L/2
    plt.figure()
    plt.plot(dop, np.array(my_depletion.thickness)/2, label='Half thickness')
    plt.plot(dop, my_depletion.depl_length, label='Depletion length')
    plt.ylabel('Depletion length, nm')
    plt.xlabel('Dopant molar fraction')
    plt.xscale('log')
    plt.legend()
    plt.savefig("Depletion_length_vs_Thickness.png")
    plt.close()

    # Fig. 6: coulomb vs. x
    plt.figure(figsize=(6, 3))
    for i_d in range(0, my_depletion.n_d):
        plt.plot(np.array(my_depletion.x1[i_d])/my_depletion.thickness[i_d], np.array(my_depletion.coulomb1[i_d]), label='doping = {:2.2E}'.format(my_depletion.dop[i_d]))
    plt.ylabel('Coulomb potential, eV')
    plt.xlabel('x, nm')
    plt.legend()
    plt.savefig("Coulomn_vs_x.png")
    plt.close()

    """
    N-1     PRINT
    """

    # uncomment to print

    print("Depletion: \n {} \n".format(my_depletion.depl_length))
    print("Mobile fraction: \n {} \n".format(my_depletion.mobile_fraction))
    print("Thickness: \n {} \n".format(my_depletion.thickness))

    """
    N       PLOT
    """


class Depletion:
    def __init__(self, n_d, n_r):
        self.n_d = n_d
        self.n_r = n_r
        self.depl_length = np.zeros(n_d).tolist()
        self.mobile_fraction = np.zeros(n_d).tolist()
        self.std_err_de = np.zeros(n_d).tolist()
        self.std_err_mo = np.zeros(n_d).tolist()
        self.thickness = np.zeros(n_d).tolist()
        self.dop = np.zeros(n_d).tolist()
        self.x1 = []     # list
        self.coulomb1 = []   # list

    def compute(self):
        fn = "dop_{}/r_{}/output_job_0"  # ugly
        depl_length, mobile_fraction = np.zeros(self.n_r).tolist(), np.zeros(self.n_r).tolist()
        for i_d in range(0, self.n_d):
            for i_r in range(0, self.n_r):
                print("read-in i_dop = {}, i_r = {}".format(i_d, i_r))
                depl_length[i_r] = extract.extract_float(fn.format(i_d, i_r), "depletion length: ")
                mobile_fraction[i_r] = extract.extract_float(fn.format(i_d, i_r), "mobile carrier fraction: ")
            #self.depl_length[i_d] = gmean(-depl_length)
            #self.mobile_fraction[i_d] = gmean(mobile_fraction)

            self.depl_length[i_d] = np.mean(-np.array(depl_length)).tolist()
            self.mobile_fraction[i_d] = np.mean(np.array(mobile_fraction)).tolist()
            self.std_err_de[i_d] = (np.std(np.array(depl_length)) / self.n_r).tolist()
            self.std_err_mo[i_d] = (np.std(np.array(mobile_fraction)) / self.n_r).tolist()
            self.dop = np.loadtxt('doping.txt').tolist()  # TODO read from the yaml of r_0

    def extract_yaml(self):
        """
        extract the thickness and other parameters from ini yaml
        :return: add attribute to, save the file
        """
        name = 'dop_{}/r_0/dl.yml'
        for i_d in range(0, self.n_d):
            fid = open(name.format(i_d), "r")
            print("read yaml for dop = {}".format(i_d))
            yaml1 = YAML(typ="safe", pure=True)
            settings = yaml1.load(fid)
            fid.close()
            self.thickness[i_d] = settings['layers'][0]['thickness']
            print('thickness = {}'.format(self.thickness[i_d]))

        name = 'dop_{}/r_{}/experiments/coulomb_average_0.dat'
        # self.x = {}     # dictionary
        # self.coulomb = {}   # dictionary
        #
        # for i_d in range(0, self.n_d):
        #     tmp = np.loadtxt(name.format(i_d, 0))[0]
        #     self.x[i_d] = tmp.tolist()  # coord
        #     self.coulomb[i_d] = np.loadtxt(name.format(i_d, 0))[1]  # C pot along
        #     for i_r in range(1, self.n_r):
        #         self.coulomb[i_d] += np.loadtxt(name.format(i_d, i_r))[1]  # C pot along
        #         print('I read coulomb | dop {} r {}'.format(i_d, i_r))
        #     self.coulomb[i_d] /= c.elementary_charge*self.n_r  # J --> eV
        # for i_d in range(0, self.n_d):
        #     self.coulomb[i_d] = self.coulomb[i_d].tolist()


        for i_d in range(0, self.n_d):
            tmp = np.loadtxt(name.format(i_d, 0))[0].tolist()
            self.x1.append(tmp)  # coord
            self.coulomb1.append(np.loadtxt(name.format(i_d, 0))[1])  # C pot along
            for i_r in range(1, self.n_r):
                self.coulomb1[i_d][:] += np.loadtxt(name.format(i_d, i_r))[1]  # C pot along
                print('I read coulomb1 | dop {} r {}'.format(i_d, i_r))
            self.coulomb1[i_d][:] /= c.elementary_charge/self.n_r  # J --> eV
        for i_d in range(0, self.n_d):
            self.coulomb1[i_d] = self.coulomb1[i_d].tolist()

    def extract_level(self, level='IP_x', mol_num='0'):

        for i_d in range(0, self.n_d):
            for i_r in range(0, self.n_r):
                name = 'dop_{}/r_{}/experiments/tmp/{}_{}_0.dat'.format(i_d, i_r, level, mol_num)
                print("I make doping {} replica {}".format(i_d, i_r))
                print("dir name = {}".format(name))
                IP_x = np.loadtxt(name)
        pass




def plot_trend(x, y, color_num = '0'):
    z = np.polyfit(x, y, 1)
    p = np.poly1d(z)
    plt.plot(x, p(x), color='C{}'.format(color_num), linestyle=':')




if __name__ == '__main__':
    main()