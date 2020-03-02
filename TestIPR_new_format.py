"""
Test IP(R) as implemented by Timo
"""

from __future__ import print_function
from __future__ import absolute_import
import numpy as np
import matplotlib
matplotlib.use('pdf')
import matplotlib.pyplot as plt
import os
import re
import yaml
import scipy.constants as const
import time

def main():

    my_mol_name = '7b043a1428ff586ba119947b48219969' #TPD

    my_mol_name = 'f273d730b77306e88fc77abb9fa56671' #aNPD

    my_qp_output = QPOutput(my_mol_name)

    folders, radii, dict_radii_folder = my_qp_output.folders, my_qp_output.radii, my_qp_output.dict_radii_folder

    # IP of R single delta mean

    my_qp_output.extract_mean_energies(mol_name = my_mol_name)

    my_qp_output.plot_single_delta()
    my_qp_output.plot_full_env()
    my_qp_output.plot_single_delta()
    my_qp_output.extract_eps()
    my_qp_output.save()



    print(my_qp_output.mean_full_env)
    print(my_qp_output.radii)


    print("I am done")

    pass

############ FUNCTIONS ################
class QPOutput:
    def __init__(self, my_mol_name):
        self.return_target_folders(my_mol_name)

    def return_target_folders(self, my_mol_name):
        folder_list = os.listdir("Analysis")

        _pattern = 'IP_'
        r = re.compile(_pattern + my_mol_name + '_for_radius_' '[0-9]*\.[0-9]*' + '_classical_correction')

        match_folders = [folder for folder in folder_list if r.match(folder)]
        #r = re.compile('.*'+'_classical_correction')
        #match_folders = [folder for folder in match_folders if r.match(folder)]

        radius_pattern = re.compile("\d+\.\d+")
        radii = np.empty(len(match_folders))
        my_dict = {}
        for i, folder in enumerate(match_folders):
            radii[i] = float(re.findall(radius_pattern, folder)[0])
            my_dict[radii[i]] = folder
        radii = np.sort(radii)

        self.folders = match_folders
        self.radii = radii
        self.dict_radii_folder = my_dict

    def extract_mean_energies(self, mol_name):
        self.mean_single_delta = np.zeros(len(self.folders))
        self.mean_full_env = np.zeros(len(self.folders))

        for i, radius in enumerate(self.radii):
            path2file = 'Analysis/' + self.dict_radii_folder[radius] + '/IP_' + mol_name + '_summary.yml'
            fid = open(path2file)
            this_dict = yaml.load(fid, Loader=yaml.SafeLoader)
            fid.close()
            self.mean_single_delta[i] = this_dict['mean_single_delta']
            self.mean_full_env[i] = this_dict['mean_full_env']
    def extract_eps(self, ab = [0, 50]):
        # fit the slope --> get the eps
        a, b = ab[0], ab[1]
        C = 7.1997176734999995  # coefficient in the formula

        target_i = np.where(np.logical_and(self.radii > a, self.radii < b))[0]
        all_r = self.radii[target_i]
        my_all_ips = -self.mean_single_delta[target_i]

        print('ips', my_all_ips)
        coef_poly = np.polyfit(1.0 / all_r, my_all_ips, 1)
        print("coef_poly = ", coef_poly)
        print("Extracted dielectric permittivity:", C / (C - coef_poly[0]))

        plt.plot(10/self.radii, -self.mean_single_delta, LineStyle='-', marker='o')
        plt.plot(10/all_r, coef_poly[0]*1/all_r + coef_poly[1])
        plt.xlabel('10/R, A')
        plt.ylabel('IP, eV')
        plt.savefig('IP_vs_1_over_R_sd_epsilon.png')
        plt.close()

    def plot_single_delta(self):
        plt.plot(self.radii, -self.mean_single_delta, LineStyle='-', marker='o')
        plt.xlabel('R, A')
        plt.ylabel('IP, eV')
        plt.savefig('IP_vs_R_sd.png')
        plt.close()

        plt.plot(10/self.radii, -self.mean_single_delta, LineStyle='-', marker='o')
        plt.xlabel('10/R, A')
        plt.ylabel('IP, eV')
        plt.savefig('IP_vs_1_over_R_sd.png')
        plt.close()

    def plot_full_env(self):

        plt.plot(self.radii, -self.mean_full_env, LineStyle='-', marker='o')
        plt.xlabel('R, A')
        plt.ylabel('IP, eV')
        plt.savefig('IP_vs_R_fe.png')
        plt.close()

        plt.plot(10/self.radii, -self.mean_full_env, LineStyle='-', marker='o')
        plt.xlabel('10/R, A')
        plt.ylabel('IP, eV')
        plt.savefig('IP_vs_1_over_R_fe.png')
        plt.close()

    def save(self):

        print("I save data to IP.dat and R.dat")
        np.savetxt('r.dat', self.radii)
        np.savetxt('IP.dat', -self.mean_single_delta)


if __name__ == '__main__':
    main()