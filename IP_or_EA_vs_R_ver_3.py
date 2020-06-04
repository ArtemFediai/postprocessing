"""
Compute/plot:
IP(R)
EA(R)
0.5(IP(R) - EA(R))
as extracted from QP (damped shells)
"""

from __future__ import print_function
from __future__ import absolute_import
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
        qp_settings_file = 'settings_ng.yml'  # default name of the QP settings file
        with open(qp_settings_file, 'r') as fid:
            qp_settings = yaml.load(fid, Loader=yaml.FullLoader)
        R = int(qp_settings['System']['Shells']['0']['cutoff'])
        ab = [R // 2, R]  # ab is the left/right limit of the epsilon evaluation interval
        my_mol_name = os.listdir('Analysis')[0].split('_')[1]  # will work only for 1-component material

        # here you can re-define R, ab, my_mol_name if necessary -->
        # R = 40
        # ab = [20, 40]
        # my_mol_name ='0ea4d81ac970d3f4fdbbe46acd91a041'
        # <--

        # get ea/ip output
        Analysis_output = QPOutput(my_mol_name)  # ini: returns target folders

        Analysis_output.extract_energies()

        Analysis_output.extract_true_vacuum_energies()

        Analysis_output.get_polarization_elementwise()

        Analysis_output.plot_p()

    print("\nI am done")
    print("Total Computation Time: {} sec".format(t.interval))


############################
class QPOutput:
    def __init__(self, mol_name):
        self.mol_name = mol_name
        self.folders_IP,\
        self.radii,\
        self.dict_radii_folder_IP = self.return_analysis_folders(IP_or_EA='IP')
        self.folders_EA, \
        self.radii,\
        self.dict_radii_folder_EA = self.return_analysis_folders(IP_or_EA='EA')
        self.inv_radii = 1.0/np.asarray(self.radii)
        self.extract_core_ids()
        self.n_mol = len(self.core_ids)

        self.E_sd_plus = {}  # E0 + V(core-env)
        self.E_sd_minus = {}
        self.E_sd_0 = {}

        self.E_sum_env_plus = {}  # E(env) + 2*V(env-env) + V(core-env)
        self.E_sum_env_minus = {}
        self.E_sum_env_0 = {}

        # 4 components of the polarization energy -->
        self.E0_plus = {}  # E0
        self.E0_minus = {} # (DFT)
        self.E0_0 = {}

        self.V0_plus = {}  # V(core-environment) < R
        self.V0_minus = {} # Coulomb
        self.V0_0 = {}

        self.E_env_plus = {}  # E(environment) < R
        self.E_env_minus = {} # DFT
        self.E_env_0 = {}

        self.V_env_plus = {}  # V(env-env) < R
        self.V_env_minus = {} # Coulomb
        self.V_env_0 = {}
        # <--


        #self.extract_energies()

        # --> 5-th component of the polarization energy
        self.E0_plus_vacuum = {}  # E0
        self.E0_minus_vacuum = {} # (DFT)
        self.E0_0_vacuum = {} # true vacuum

        self.vacuum_folder_path = './vacuum' #where to look it for
        #self.extract_true_vacuum_energies()
        # <--


    def return_analysis_folders(self, IP_or_EA):  #
        # return relevant folders in Analysis
        # return radii
        # return radii and folders as the dictionary
        folder_list = os.listdir("Analysis")
        r = re.compile(IP_or_EA + '_' + self.mol_name + '_for_radius_' '[0-9]*\.[0-9]*' + '_classical_correction')
        match_folders = [folder for folder in folder_list if r.match(folder)]
        radius_pattern = re.compile("\d+\.\d+")
        radii = np.empty(len(match_folders))
        my_dict = {}
        for i, folder in enumerate(match_folders):
            radii[i] = float(re.findall(radius_pattern, folder)[0])
            my_dict[radii[i]] = folder
        radii = np.sort(radii)
        return match_folders, radii, my_dict


    def extract_core_ids(self):
        """
        extracts all core ides as a list of strings
        """
        path2folder = 'Analysis/' + self.dict_radii_folder_IP[self.radii[0]]
        dirs = [dir for dir in os.listdir(path2folder) if dir.endswith('1.yml')]
        self.core_ids = [int(this_dir.split('_')[2]) for this_dir in dirs]

    def extract_energies(self):
        """
        extract all components of the polarization (except for the true vacuum)
        """


        for i, radius in enumerate(self.radii):
            self.E_sd_plus[radius] = {}
            self.E_sd_0[radius] = {}
            self.E_sd_minus[radius] = {}

            self.E_sum_env_plus[radius] = {}
            self.E_sum_env_0[radius] = {}
            self.E_sum_env_minus[radius] = {}

            self.V0_plus[radius] = {}
            self.V0_0[radius] = {}
            self.V0_minus[radius] = {}

            self.E_env_plus[radius] = {}
            self.E_env_0[radius] = {}
            self.E_env_minus[radius] = {}

            self.V_env_plus[radius] = {}
            self.V_env_0[radius] = {}
            self.V_env_minus[radius] = {}

            for j, core_id in enumerate(self.core_ids):
                path2file_ip = \
                    'Analysis/' + self.dict_radii_folder_IP[radius] + '/Matrix-analysis-IP_' \
                    + self.mol_name + '-Mol_' + str(core_id) + '_C_1.yml'
                path2file_ea = \
                    'Analysis/' + self.dict_radii_folder_EA[radius] + '/Matrix-analysis-EA_' \
                    + self.mol_name + '-Mol_' + str(core_id) + '_C_-1.yml'

                # IP. Charged states: "+" and "0"
                with open(path2file_ip) as fid:
                    ip_dict = yaml.load(fid, Loader=yaml.SafeLoader)
                with open(path2file_ea) as fid:
                    ea_dict = yaml.load(fid, Loader=yaml.SafeLoader)


                # sd extraction. E_sd = E_0 + V_0
                self.E_sd_plus[radius][core_id] = ip_dict[int(core_id)]['total_e_charged']
                self.E_sd_0[radius][core_id] = ip_dict[core_id]['total_e_uncharged']
                self.E_sd_minus[radius][core_id] = ea_dict[int(core_id)]['total_e_charged']
                # E_0
                self.E0_plus[core_id] = ip_dict[int(core_id)]['total_e_charged_vacuum']
                self.E0_0[core_id] = ip_dict[int(core_id)]['total_e_uncharged_vacuum']
                self.E0_minus[core_id] = ea_dict[int(core_id)]['total_e_charged_vacuum']
                # V_0
                self.V0_plus[radius][core_id] = self.E_sd_plus[radius][core_id] - self.E0_plus[core_id]
                self.V0_0[radius][core_id] = self.E_sd_0[radius][core_id] - self.E0_0[core_id]
                self.V0_minus[radius][core_id] = self.E_sd_minus[radius][core_id] - self.E0_minus[core_id]

                # E_sum_env = \sum_i\ne 0  E_i  \sum_{j=0}^{N} V_{ij}
                ip_env_sub_dict = ip_dict[int(core_id)]['environment_molecules']
                ea_env_sub_dict = ea_dict[int(core_id)]['environment_molecules']

                list_total_e_env_plus = [ip_env_sub_dict[env_id]['total_e_charged'] for env_id in ip_env_sub_dict]
                self.E_sum_env_plus[radius][int(core_id)] = np.sum(list_total_e_env_plus) if not list_total_e_env_plus == [] else 0.0
                list_total_e_env_0 = [ip_env_sub_dict[env_id]['total_e_uncharged'] for env_id in ip_env_sub_dict]
                self.E_sum_env_0[radius][int(core_id)] = np.sum(list_total_e_env_0) if not list_total_e_env_0 == [] else 0.0
                list_total_e_env_minus = [ea_env_sub_dict[env_id]['total_e_charged'] for env_id in ea_env_sub_dict]
                self.E_sum_env_minus[radius][int(core_id)] = np.sum(list_total_e_env_minus) if not list_total_e_env_minus == [] else 0.0

                # E_env = \sum_i \ne 0  E_i. sum of DFT env energies.
                list_vacuum_env_e_plus = [ip_env_sub_dict[env_id]['total_e_charged_vacuum'] for env_id in ip_env_sub_dict]
                self.E_env_plus[radius][int(core_id)] = np.sum(list_vacuum_env_e_plus) if not list_vacuum_env_e_plus == [] else 0.0
                list_vacuum_env_e_0 = [ip_env_sub_dict[env_id]['total_e_uncharged_vacuum'] for env_id in ip_env_sub_dict]
                self.E_env_0[radius][int(core_id)] = np.sum(list_vacuum_env_e_0) if not list_vacuum_env_e_0 == [] else 0.0
                list_vacuum_env_e_minus = [ea_env_sub_dict[env_id]['total_e_charged_vacuum'] for env_id in ea_env_sub_dict]
                self.E_env_minus[radius][int(core_id)] = np.sum(list_vacuum_env_e_minus) if not list_vacuum_env_e_minus == [] else 0.0

                # V_env = 0.5 (\sum_{i=1} \sum_{j=1} V_{ij}). classical interaction of env. mols
                self.V_env_plus[radius][core_id] = 0.5 * (self.E_sum_env_plus[radius][core_id]
                                                          - self.E_env_plus[radius][core_id]
                                                          - self.V0_plus[radius][core_id])

                self.V_env_0[radius][core_id] = 0.5 * (self.E_sum_env_0[radius][core_id]
                                                          - self.E_env_0[radius][core_id]
                                                          - self.V0_0[radius][core_id])

                self.V_env_minus[radius][core_id] = 0.5 * (self.E_sum_env_minus[radius][core_id]
                                                          - self.E_env_minus[radius][core_id]
                                                          - self.V0_minus[radius][core_id])

        with open('Analysis/energies.pkl', 'wb') as fid:
            pickle.dump([self.E0_plus, self.E0_0, self.E0_minus,
                         self.V0_plus, self.V0_0, self.V0_minus,
                         self.V_env_plus, self.V_env_0, self.V_env_minus,
                         self.E_env_plus, self.E_env_0, self.E_env_minus],
                         fid)
        # TODO delete to check if exists above and do load without doing
        with open('Analysis/energies.pkl', 'rb') as fid:
            [self.E0_plus, self.E0_0, self.E0_minus,
             self.V0_plus, self.V0_0, self.V0_minus,
             self.V_env_plus, self.V_env_0, self.V_env_minus,
             self.E_env_plus, self.E_env_0, self.E_env_minus] \
                = pickle.load(fid)


        append_dict_with_mean(self.V0_plus, self.V0_0, self.V0_minus,
                              self.V_env_plus, self.V_env_0, self.V_env_minus,
                              self.E_env_plus, self.E_env_0, self.E_env_minus,
                              self.E0_plus, self.E0_0, self.E0_minus)  # compute and add "mean" to all mentioned dicts

        print('here')

    def extract_true_vacuum_energies(self):
        if os.path.exists(self.vacuum_folder_path):
            ip_vacuum_dir = self.vacuum_folder_path + '/Analysis/' + self.folders_IP[0].split('_')[0]  + '_' + self.folders_IP[0].split('_')[1]
            ea_vacuum_dir = self.vacuum_folder_path + '/Analysis/' + self.folders_EA[0].split('_')[0]  + '_' + self.folders_EA[0].split('_')[1]
            for core_id in self.core_ids:
                with open(ip_vacuum_dir+'/Matrix-analysis-IP_'+self.mol_name+ '-Mol_'+str(core_id) + '_C_1.yml') as fid:
                    ip_dict = yaml.load(fid, Loader=yaml.SafeLoader)
                with open(ea_vacuum_dir+'/Matrix-analysis-EA_'+self.mol_name+ '-Mol_'+str(core_id) + '_C_-1.yml') as fid:
                    ea_dict = yaml.load(fid, Loader=yaml.SafeLoader)

                self.E0_plus_vacuum[core_id] = ip_dict[core_id]['total_e_charged_vacuum']
                self.E0_0_vacuum[core_id] = ip_dict[core_id]['total_e_uncharged_vacuum']
                self.E0_minus_vacuum[core_id] = ea_dict[core_id]['total_e_charged_vacuum']
            append_dict_with_mean(self.E0_plus_vacuum, self.E0_0_vacuum, self.E0_minus_vacuum)

        else:
            self.E0_plus_vacuum = self.E0_plus
            self.E0_0_vacuum = self.E0_0
            self.E0_minus_vacuum = self.E0_minus

    def get_polarization_elementwise(self):
        e0p_v = self.E0_plus_vacuum['mean']
        e0m_v = self.E0_minus_vacuum['mean']
        e00_v = self.E0_0_vacuum['mean']

        e0p = self.E0_plus['mean']
        e0m = self.E0_minus['mean']
        e00 = self.E0_0['mean']

        v0p = [self.V0_plus[radius]['mean'] for radius in self.V0_plus]
        v0m = [self.V0_minus[radius]['mean'] for radius in self.V0_minus]
        v00 = [self.V0_0[radius]['mean'] for radius in self.V0_0]

        vep = [self.V_env_plus[radius]['mean'] for radius in self.V_env_plus]
        vem = [self.V_env_minus[radius]['mean'] for radius in self.V_env_minus]
        ve0 = [self.V_env_0[radius]['mean'] for radius in self.V_env_0]

        eep = [self.E_env_plus[radius]['mean'] for radius in self.E_env_plus]
        eem = [self.E_env_minus[radius]['mean'] for radius in self.E_env_minus]
        ee0 = [self.E_env_0[radius]['mean'] for radius in self.E_env_0]

        v0p, v0m, v00, vep, vem, ve0, eep, eem, ee0 = list2arr(v0p, v0m, v00, vep, vem, ve0, eep, eem, ee0)

        # cation
        p1 = (e0p_v - e00_v) - (e0p - e00)
        p2 = -(v0p - v00)
        p3 = -(vep - ve0)
        p4 = -(eep - ee0)

        self.dE0_cation = p1  # vacuum; core
        self.dV0_cation = p2  # core-env
        self.dVe_cation = p3  # env-env
        self.dEe_cation = p4  # env DFT
        self.P_cation = p1+p2+p3+p4

        # anion
        p1 = (e0m_v - e00_v) - (e0m - e00)
        p2 = -(v0m - v00)
        p3 = -(vem - ve0)
        p4 = -(eem - ee0)

        self.dE0_anion = p1  # vacuum; core
        self.dV0_anion = p2  # core-env
        self.dVe_anion = p3  # env-env
        self.dEe_anion = p4  # env DFT

        self.P_anion = p1+p2+p3+p4


    def plot_p(self, xlim = [0, 1.25]):
        plt.plot(10*self.inv_radii, self.dE0_cation*np.ones(len(self.radii)), label = '$ \Delta E_{core}$ ')
        plt.plot(10*self.inv_radii, self.dV0_cation, label = '$\Delta V_{core}$')
        plt.plot(10*self.inv_radii, self.dVe_cation, label = '$\Delta V_{env}$')
        plt.plot(10*self.inv_radii, self.dEe_cation, label = '$\Delta E_{env}$')
        plt.plot(10*self.inv_radii, self.P_cation, label = '$P$')
        plt.xlabel('$10/R,  10 / \AA^{-1} $')
        plt.ylabel('$P^+$, eV')
        plt.legend()
        plt.xlim(xlim)
        plt.grid()
        ax1 = plt.gca()
        add_inverse_axis(ax1)
        plt.savefig('P.png')


def add_inverse_axis(initial_axis, rs_plot=np.array([1, 2, 3, 4, 5, 7, 10, 20, 30, 40, 50]),
                     rs_grid=np.array([])):
    def inv(x):
        return 10 / x

    ax2 = initial_axis.twiny()
    ax2.set_xlabel('$R, \AA$')
    ax2.set_xticks(inv(rs_plot), minor=False)
    ax2.set_xticks(inv(rs_grid), minor=True)
    ax2.set_xticklabels(rs_plot)
    ax2.set_xbound(initial_axis.get_xbound())
    ax2.grid(True, which='minor', color = 'blue')
    ax2.grid(True, which='major', color = 'red', linestyle='dotted')


class Timer:
    def __enter__(self):
        self.start = time.process_time()
        return self

    def __exit__(self, *args):
        self.end = time.process_time()
        self.interval = self.end - self.start


def append_folder(paths_to_mols, max_or_min):
    """
    append the list of strings (which represent paths) with the last first number in these paths.
    numbers are nested folders in path.
    :param paths_to_mols: array of strings that mean paths to QP output for molecules
    :param max_or_min: last or first step
    :return: updated paths to mols
    """

    paths_to_mols = np.atleast_1d(paths_to_mols)
    all_steps = [int(number) for number in
                 os.listdir(paths_to_mols[0])]  # there must be the same number of steps for every molecule
    if max_or_min == 'max':
        step = str(max(all_steps))
    else:
        step = str(min(all_steps))
    paths_to_mols = [path + '/' + step for path in paths_to_mols]
    return paths_to_mols

def append_dict_with_mean(*dictionaries):
    """
    d = [1: {1:2, 2:3}, 2: {1:2, 2:3}]
    :param d:
    :return:
    d = [1: {1:2, 2:3, 'mean': 2.5}, 2: {1:2, 2:3, 'mean': 2.5}]
    """


    for d in dictionaries:
        my_depth = depth(d)
        if my_depth > 1:
            raise UserWarning("I cannot handle dictionaries which are > 1 deep")
        # print('depth of the vocabulary: ', my_depth)
        if not isinstance(list(d.values())[0], dict):
            this_mean = np.mean(list(d.values()))
            d['mean'] = this_mean
        else:
            for key, value in d.items():
                this_mean = np.mean(list(value.values()))
                d[key]['mean'] = this_mean

def depth(d):
    depth = 0
    q = [(i, depth + 1) for i in d.values() if isinstance(i, dict)]
    max_depth = 0
    while (q):
        n, depth = q.pop()
        max_depth = max(max_depth, depth)
        # q = q + [(i, depth + 1) for i in n.values() if isinstance(i, dict)]
    return max_depth

def list2arr(*lists):

    new_list = []
    for list in lists:
        new_list.append(np.asarray(list))
    output = tuple(new_list)
    return output

if __name__ == '__main__':
    main()
