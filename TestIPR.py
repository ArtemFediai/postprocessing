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

    my_mol_name = '7b043a1428ff586ba119947b48219969'


    my_qp_output = QPOutput()

    folders, radii, dict_radii_folder = my_qp_output.folders, my_qp_output.radii, my_qp_output.dict_radii_folder

    # IP of R single delta mean

    my_qp_output.extract_mean_energies(mol_name = my_mol_name)

    my_qp_output.plot_single_delta()
    my_qp_output.plot_full_env()
    my_qp_output.plot_single_delta()



    print(my_qp_output.mean_full_env)
    print(my_qp_output.radii)


    print("I am done")

    pass

############ FUNCTIONS ################
class QPOutput:
    def __init__(self):
        self.return_target_folders()

    def return_target_folders(self):
        _pattern = 'EAIP_for_radius_'
        r = re.compile(_pattern + '[0-9]*\.[0-9]*')

        folder_list = os.listdir("Analysis")
        match_folders = [folder for folder in folder_list if r.match(folder)]
        r = re.compile('.*'+'_classical_correction')
        match_folders = [folder for folder in match_folders if not r.match(folder)]

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
            path2file = 'Analysis/' + self.dict_radii_folder[radius] + '/eaip_' + mol_name + '_IP_summary.yml'
            fid = open(path2file)
            this_dict = yaml.load(fid, Loader=yaml.SafeLoader)
            fid.close()
            self.mean_single_delta[i] = this_dict['mean_single_delta']
            self.mean_full_env[i] = this_dict['mean_full_env']

    def plot_single_delta(self):
        plt.plot(self.radii, -self.mean_single_delta, LineStyle='-', marker='o')
        plt.xlabel('R, A')
        plt.ylabel('IP, eV')
        plt.savefig('IP_vs_R_sd.png')
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




if __name__ == '__main__':
    main()