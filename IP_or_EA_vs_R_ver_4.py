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
        R = int(qp_settings['System']['Shells']['0']['cutoff']) # this will work only for one shell!
        ab = [19, 31]  # ab is the left/right limit of the epsilon evaluation interval
        number_density_nm3_aNPD = 1.158623303810171 # aNPD
        number_density_nm3_TCTA = 0.9162165860886232 # TCTA
        number_density_nm3_C60 = 1.4266923 # C60
        number_density_nm3_SpiroOMeTAD = 0.4873207140786595  # SpiroOMeTAD
        #
        my_density = number_density_nm3_SpiroOMeTAD
        #
#        my_mol_name = os.listdir('Analysis')[0].split('_')[1]  # will work only for 1-component material #TODO: wrong

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

        Analysis_output.plot_p_n()

        Analysis_output.plot_IPEA()

        Analysis_output.extract_eps_from_polarization(ab=ab, number_density_nm3=my_density)

        Analysis_output.extract_eps_from_polarization_new(ab=ab, number_density_nm3=my_density)

        Analysis_output.save_polarization_elementwise()

    print("\nI am done")
    print("Total Computation Time: {} sec".format(t.interval))


############################
class QPOutput:
    def __init__(self):
        self.mol_name = []
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

        # --> n_mols
#        self.n_mols = np.zeros(len(self.radii), dtype=np.int)
        self.n_mols = {}
        # <-- n_mols

        # mean polarization energies that will computed and saved
        self.dE0_anion = []  # vacuum; core
        self.dV0_anion = []  # core-env
        self.dVe_anion = []  # env-env
        self.dEe_anion = []  # env DFT
        self.P_anion = []

        self.dE0_cation = []  # vacuum; core
        self.dV0_cation = []  # core-env
        self.dVe_cation = []  # env-env
        self.dEe_cation = []  # env DFT
        self.P_cation = []

        self.radii_renormalized = []

    def return_analysis_folders(self, IP_or_EA):  #
        # return relevant folders in Analysis
        # return radii
        # return radii and folders as the dictionary
        folder_list = os.listdir("Analysis/" + IP_or_EA + '_by_radius')
        r = re.compile(IP_or_EA + '_for_radius_' '[0-9]*\.[0-9]*' + '_classical_correction')
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
        path2folder = 'Analysis/IP_by_radius/' + self.dict_radii_folder_IP[self.radii[0]] + '/'
        analysis_files = [dir for dir in os.listdir(path2folder) if dir.startswith('Matrix-analysis-IP_')]
        analysis_file = path2folder + analysis_files[0]  #work for 1 component system
        with open(analysis_file, 'r') as fid:
            my_file = yaml.load(fid, Loader=yaml.FullLoader)
        self.core_ids = list(my_file.keys())
        self.mol_name = analysis_files[0].split('_')[1].split('.')[0]


        print('coreids', self.core_ids)

    def extract_energies(self):
        """
        extract all components of the polarization (except for the true vacuum)
        """
        path2save = 'Analysis/energies.pkl'
        #check, if I have to extract them, or they are already extracted. This the latter case, load them.
        if os.path.exists(path2save):
            print("extraction of the polarizaion has already been done. Loading polarizations from from pkl")
            # TODO delete to check if exists above and do load without doing
            with open('Analysis/energies.pkl', 'rb') as fid:
                [self.E0_plus, self.E0_0, self.E0_minus,
                 self.V0_plus, self.V0_0, self.V0_minus,
                 self.V_env_plus, self.V_env_0, self.V_env_minus,
                 self.E_env_plus, self.E_env_0, self.E_env_minus,
                 self.n_mols] \
                    = pickle.load(fid)
        else:
            print('Energies are being extracting and will be saved to pkl')
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

                self.n_mols[radius] = {}

                for j, core_id in enumerate(self.core_ids):
                    #path2file_ip = \
                    #    'Analysis/' + self.dict_radii_folder_IP[radius] + '/Matrix-analysis-IP_' \
                    #    + self.mol_name + '-Mol_' + str(core_id) + '_C_1.yml'

                    path2file_ip = \
                        'Analysis/IP_by_radius/' + self.dict_radii_folder_IP[radius]\
                        + '/Matrix-analysis-IP_' + self.mol_name + '.yml' # new
                    path2file_ea = \
                        'Analysis/EA_by_radius/' + self.dict_radii_folder_EA[radius]\
                        + '/Matrix-analysis-EA_' + self.mol_name + '.yml'

                    # IP. Charged states: "+" and "0"
                    with open(path2file_ip) as fid:
                        ip_dict = yaml.load(fid, Loader=yaml.SafeLoader)
                    with open(path2file_ea) as fid:
                        ea_dict = yaml.load(fid, Loader=yaml.SafeLoader)


                    # number of mols extraction
                    self.n_mols[radius][core_id] = len(ip_dict[int(core_id)]['energies'])

                    # sd extraction. E_sd = E_0 + V_0
                    self.E_sd_plus[radius][core_id] = ip_dict[int(core_id)]['energies'][int(core_id)]['total_e_charged']  #new
                    self.E_sd_0[radius][core_id] = ip_dict[core_id]['energies'][int(core_id)]['total_e_uncharged']
                    self.E_sd_minus[radius][core_id] = ea_dict[int(core_id)]['energies'][int(core_id)]['total_e_charged']
                    # E_0
                    self.E0_plus[core_id] = ip_dict[int(core_id)]['energies'][int(core_id)]['total_e_charged_vacuum']
                    self.E0_0[core_id] = ip_dict[int(core_id)]['energies'][int(core_id)]['total_e_uncharged_vacuum']
                    self.E0_minus[core_id] = ea_dict[int(core_id)]['energies'][int(core_id)]['total_e_charged_vacuum']
                    # # E_0_vacuum
                    # self.E0_plus_vacuum[core_id] =
                    # self.E0_0_vacuum[core_id] =
                    # self.E0_minus_vacuum[core_id] =


                    # V_0
                    self.V0_plus[radius][core_id] = self.E_sd_plus[radius][core_id] - self.E0_plus[core_id]
                    self.V0_0[radius][core_id] = self.E_sd_0[radius][core_id] - self.E0_0[core_id]
                    self.V0_minus[radius][core_id] = self.E_sd_minus[radius][core_id] - self.E0_minus[core_id]

                    # E_sum_env = \sum_i\ne 0  E_i  \sum_{j=0}^{N} V_{ij}
                    ip_env_sub_dict = ip_dict[int(core_id)]['energies']#new
                    del ip_env_sub_dict[int(core_id)]
    #                del ip_env_sub_dict['info']  # TODO: do I need to dlt this?


                    ea_env_sub_dict = ea_dict[int(core_id)]['energies']  # new
                    del ea_env_sub_dict[int(core_id)]
    #                del ea_env_sub_dict['info']  # TODO: do I need to dlt this?

    #                tmp = ip_env_sub_dict['energies'][]

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


            append_dict_with_mean(self.V0_plus, self.V0_0, self.V0_minus,
                                  self.V_env_plus, self.V_env_0, self.V_env_minus,
                                  self.E_env_plus, self.E_env_0, self.E_env_minus,
                                  self.E0_plus, self.E0_0, self.E0_minus,
                                  self.n_mols)  # compute and add "mean" to all mentioned dicts

            with open('Analysis/energies.pkl', 'wb') as fid:
                pickle.dump([self.E0_plus, self.E0_0, self.E0_minus,
                             self.V0_plus, self.V0_0, self.V0_minus,
                             self.V_env_plus, self.V_env_0, self.V_env_minus,
                             self.E_env_plus, self.E_env_0, self.E_env_minus,
                             self.n_mols],
                             fid)
            print("Energies are extracted and dumped to pkl")





    def extract_true_vacuum_energies(self):
        if os.path.exists('Analysis/vacuum_energies.pkl'):
            with open('Analysis/vacuum_energies.pkl', 'rb') as fid:
                print("extraction (vacuum) has already been done. Loading p.e. from pkl")
                [self.E0_plus_vacuum, self.E0_0_vacuum, self.E0_minus_vacuum] = pickle.load(fid)
        else:
            print('Vacuum energies are being extracting and will be saved to pkl')
            zero_dir = 'quantumpatch_runtime_files/uncharged'
            zero_file = zero_dir + '/' + str((min([int(this_dir) for this_dir in os.listdir(zero_dir)]))) + '/energies.ene.yml.gz'
            zero_dict = load_yaml_from_gzip(zero_file)
            for core_id in self.core_ids:
                plus_dir = 'quantumpatch_runtime_files/Mol_' + str(core_id) + '_C_1'
                minus_dir = 'quantumpatch_runtime_files/Mol_' + str(core_id) + '_C_-1'
                plus_file = plus_dir + '/' + str((min([int(this_dir) for this_dir in os.listdir(plus_dir)]))) + '/energies.ene.yml.gz'
                minus_file = minus_dir + '/' + str((min([int(this_dir) for this_dir in os.listdir(minus_dir)])))  + '/energies.ene.yml.gz'
                minus_dict = load_yaml_from_gzip(minus_file)
                plus_dict = load_yaml_from_gzip(plus_file)

                self.E0_plus_vacuum[core_id] = plus_dict[str(core_id)]['total']
                self.E0_0_vacuum[core_id] = zero_dict[str(core_id)]['total']
                self.E0_minus_vacuum[core_id] = minus_dict[str(core_id)]['total']

            append_dict_with_mean(self.E0_plus_vacuum, self.E0_0_vacuum, self.E0_minus_vacuum)
            with open('Analysis/vacuum_energies.pkl', 'wb') as fid:
                pickle.dump([self.E0_plus_vacuum, self.E0_0_vacuum, self.E0_minus_vacuum], fid)

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

    def save_polarization_elementwise(self):
        f4p = 'Analysis/p_elementwise/'
        if not os.path.exists(f4p):
            os.mkdir(f4p)

        np.savetxt(fname=f4p+'dE0_cation.dat',  X=self.dE0_cation*np.ones(len(self.radii)))
        np.savetxt(fname=f4p+'dE0_anion.dat',  X=self.dE0_anion*np.ones(len(self.radii)))
        np.savetxt(fname=f4p+'dV0_cation.dat',  X=self.dV0_cation)
        np.savetxt(fname=f4p+'dV0_anion.dat',  X=self.dV0_anion)
        np.savetxt(fname=f4p+'dVe_cation.dat',  X=self.dVe_cation)
        np.savetxt(fname=f4p+'dVe_anion.dat',  X=self.dVe_anion)
        np.savetxt(fname=f4p+'dEe_cation.dat',  X=self.dEe_cation)
        np.savetxt(fname=f4p+'dEe_anion.dat',  X=self.dEe_anion)
        np.savetxt(fname=f4p+'radii.dat',  X=self.radii)
        np.savetxt(fname=f4p+'radii_renormalized.dat',  X=self.radii_renormalized)


    def plot_p(self, xlim = [0, 1.25]):
        # plots cation and anion polarization: total end element-wise

        x = 10*self.inv_radii

        # cation
        #Pi = [self.dE0_cation*np.ones(len(self.radii)), self.dV0_cation, self.dVe_cation, self.dEe_cation]
        Pi = [self.dE0_cation*np.ones(len(self.radii)), self.dV0_cation, self.dVe_cation, self.dEe_cation, self.dVe_cation + self.dEe_cation]
        p = self.P_cation
        #labels_Pi = ['$\Delta E_{core}$', '$\Delta V_{core}$', '$\Delta V_{env}$', '$\Delta E_{env}$']
        labels_Pi = ['$\Delta E_{core}$', '$\Delta V_{core}$', '$\Delta V_{env}$', '$\Delta E_{env}$', '$\Delta E_{env} + \Delta V_{env}$']
        label_p =  '$P$ '
        y_name = '$P^+$'
        filename = 'P_plus.png' # without extention
        plot_4_plus_1(x, Pi,p, y_name, filename, labels_Pi, label_p, xlim)

        # anion
        #Pi = [self.dE0_anion*np.ones(len(self.radii)), self.dV0_anion, self.dVe_anion, self.dEe_anion]
        Pi = [self.dE0_anion*np.ones(len(self.radii)), self.dV0_anion, self.dVe_anion, self.dEe_anion, self.dVe_anion + self.dEe_anion]
        p = self.P_anion
        y_name = '$P^-$'
        filename = 'P_minus.png' # without extention
        plot_4_plus_1(x, Pi,p, y_name, filename, labels_Pi, label_p, xlim)

    def plot_p_n(self, xlim = [0, 1.0]):
        # AS A FUNCTION OF N**(-1/3)
        # plots cation and anion polarization: total end element-wise

        # x = 10*self.inv_radii

        mean_n = np.asarray([self.n_mols[radius]['mean'] for radius in self.n_mols])
        x = mean_n**(-1/3)

        # cation
        #Pi = [self.dE0_cation*np.ones(len(self.radii)), self.dV0_cation, self.dVe_cation, self.dEe_cation]
        Pi = [self.dE0_cation*np.ones(len(self.radii)), self.dV0_cation, self.dVe_cation, self.dEe_cation, self.dVe_cation + self.dEe_cation]
        p = self.P_cation
        #labels_Pi = ['$\Delta E_{core}$', '$\Delta V_{core}$', '$\Delta V_{env}$', '$\Delta E_{env}$']
        labels_Pi = ['$\Delta E_{core}$', '$\Delta V_{core}$', '$\Delta V_{env}$', '$\Delta E_{env}$', '$\Delta E_{env} + \Delta V_{env}$']
        label_p =  '$P$ '
        y_name = '$P^+$'
        filename = 'P_plus_vs_N.png' # without extention
        plot_4_plus_1_for_p_vs_n(x, Pi,p, y_name, filename, labels_Pi, label_p, xlim)

        # anion
        #Pi = [self.dE0_anion*np.ones(len(self.radii)), self.dV0_anion, self.dVe_anion, self.dEe_anion]
        Pi = [self.dE0_anion*np.ones(len(self.radii)), self.dV0_anion, self.dVe_anion, self.dEe_anion, self.dVe_anion + self.dEe_anion]
        p = self.P_anion
        y_name = '$P^-$'
        filename = 'P_minus_vs_N.png' # without extention
        plot_4_plus_1_for_p_vs_n(x, Pi,p, y_name, filename, labels_Pi, label_p, xlim)

    def plot_IPEA(self, xlim = [0, 1.25]):
        # plots cation and anion polarization: total end element-wise

        x = 10*self.inv_radii
        # cation
        self.IP = (self.E0_plus_vacuum['mean'] - self.E0_0_vacuum['mean'])*np.ones(len(self.P_cation)) - self.P_cation #
        plot_1(x, self.IP, 'IP.png', 'IP', xlim)
        # anion
        self.EA = (self.E0_0_vacuum['mean'] - self.E0_minus_vacuum['mean'])*np.ones(len(self.P_anion)) + self.P_anion #
        plot_1(x, self.EA, 'EA.png', 'EA', xlim)


    def extract_eps_from_polarization(self, ab=[20, 60], a1b1=[5, 1E100], number_density_nm3=1.158623303810171):
        # fit the slope --> get the eps
        self.radii = np.array(self.radii)
        mean_energies = 0.5 * (self.P_cation + self.P_anion)  # mean polarization energy

        c = 7.1997176734999995  # coefficient in the formula

        a, b, a1, b1 = ab[0], ab[1], a1b1[0], a1b1[1]  # analyze interval; plot interval

        target_i = np.where(np.logical_and(self.radii > a, self.radii < b))[0]
        selected_r = self.radii[target_i]  # selected interval to compute epsilon
        selected_energies = mean_energies[target_i]  # TODO: all energies to IP/EA

        coef_poly = np.polyfit(1.0 / selected_r, selected_energies, 1)

        selected_fitted_energies = coef_poly[0] / selected_r + coef_poly[1]

        epsilon = c / (c + coef_poly[0])

        point_a = [0, coef_poly[1]]
        point_b = [10/selected_r[-1], selected_energies[-1]]

        # fictitious eps
        s_25 = c*(1-2.5)/2.5
        s_30 = c*(1-3.0)/3.0
        s_35 = c*(1-3.5)/3.5

        point_a_25 = [0, selected_energies[-1] - s_25*1/selected_r[-1]]
        point_a_30 = [0, selected_energies[-1] - s_30*1/selected_r[-1]]
        point_a_35 = [0, selected_energies[-1] - s_35*1/selected_r[-1]]
        # fictitious eps


        self.epsilon = epsilon

        print("Slope: ", coef_poly[0])
        print("Polarization energy: ", coef_poly[1])
        print("Dielectric permittivity:", epsilon)

        plt.plot(10 / self.radii, mean_energies, LineStyle='-', marker='o')  # whole curve
        plt.plot(10 / selected_r, selected_fitted_energies)  # ??

        plt.plot([point_a[0], point_b[0]],[point_a[1],point_b[1]], label='approx')
        plt.plot([point_a_25[0], point_b[0]],[point_a_25[1],point_b[1]], label='eps 2.5', linestyle=':')
        plt.plot([point_a_30[0], point_b[0]],[point_a_30[1],point_b[1]], label='eps 3.0', linestyle=':')
        plt.plot([point_a_35[0], point_b[0]],[point_a_35[1],point_b[1]], label='eps 3.5', linestyle=':')

        plt.xlabel('10/R, A')
        plt.ylabel('polarization, eV')
        plt.xlim([10 / b1, 10 / a1])
        # plt.ylim([3.2, 3.5])  # hard-coded
        ym_ym = plt.gca().get_ylim()
        plt.plot(10.0 / np.array([20, 20]), ym_ym, label='20 A')  # TODO: make it professional
        plt.plot(10.0 / np.array([30, 30]), ym_ym, label='30 A')
        plt.plot(10.0 / np.array([40, 40]), ym_ym, label='40 A')
        plt.plot(10.0 / np.array([50, 50]), ym_ym, label='50 A')
        plt.plot(10.0 / np.array([60, 60]), ym_ym, label='60 A')
        plt.legend()
        plt.savefig('EPS_FULL_ENV_CUSTOM.png', dpi=600)
        plt.close()
#
        # below: polarization vs. n**(-1/3)

        n_mean = np.asarray([self.n_mols[radius]['mean'] for radius in self.n_mols]) #mean number of molecules
        r_renormalized = (n_mean/number_density_nm3 * 3 / 4 / np.pi * 1E3)**(1/3)  # in A
        self.radii_renormalized = r_renormalized
#        plt.plot(self.radii, r_renormalized) renormalized r vs normal r
#        plt.savefig('test.png')
        #target_i is used
        selected_r_renormalized = r_renormalized[target_i]
        # selected_energies are used

        # computing eps
        coef_poly = np.polyfit(1.0 / selected_r_renormalized, selected_energies, 1)
        selected_fitted_energies = coef_poly[0] / selected_r_renormalized + coef_poly[1]
        epsilon = c / (c + coef_poly[0])
        point_a = [0, coef_poly[1]]
        point_b = [10/selected_r_renormalized[-1], selected_energies[-1]]

        self.epsilon_renormalized = epsilon
        print("Slope renormalized: ", coef_poly[0])
        print("Polarization energy renormalized: ", coef_poly[1])
        print("Dielectric permittivity renormalized:", epsilon)


    def extract_eps_from_polarization_new(self, ab=[20, 40], a1b1=[5, 1E100], number_density_nm3 = 1.15862):
        # fit the slope --> get the eps
        self.radii = np.array(self.radii)
        mean_energies = 0.5 * (self.P_cation + self.P_anion)  # mean polarization energy

        c = 7.1997176734999995  # coefficient in the formula

        a, b, a1, b1 = ab[0], ab[1], a1b1[0], a1b1[1]  # analyze interval; plot interval

        target_i = np.where(np.logical_and(self.radii > a, self.radii < b))[0]
        selected_r = self.radii[target_i]  # selected interval to compute epsilon
        selected_energies = mean_energies[target_i]  # TODO: all energies to IP/EA

        coef_poly = np.polyfit(1.0 / selected_r, selected_energies, 1)

        selected_fitted_energies = coef_poly[0] / selected_r + coef_poly[1]

        epsilon = c / (c + coef_poly[0])

        point_a = [0, coef_poly[1]]
        point_b = [10/selected_r[-1], selected_energies[-1]]

        # fictitious eps
        s_25 = c*(1-2.5)/2.5
        s_30 = c*(1-3.0)/3.0
        s_35 = c*(1-3.5)/3.5

        point_a_25 = [0, selected_energies[-1] - s_25*1/selected_r[-1]]
        point_a_30 = [0, selected_energies[-1] - s_30*1/selected_r[-1]]
        point_a_35 = [0, selected_energies[-1] - s_35*1/selected_r[-1]]
        # fictitious eps

        self.epsilon = epsilon

        print("Slope: ", coef_poly[0])
        print("Polarization energy: ", coef_poly[1])
        print("Dielectric permittivity:", epsilon)

        plt.plot(10 / self.radii, mean_energies, LineStyle='none', marker='o', color = 'C0', label = 'polarization', mfc='none')  # whole curve
        plt.plot(10 / selected_r, selected_fitted_energies)  # ??

        plt.plot([point_a[0], point_b[0]],[point_a[1],point_b[1]], label='interpolation',color='C0', LineStyle=':')


        #

        # below: polarization vs. n**(-1/3)

        n_mean = np.asarray([self.n_mols[radius]['mean'] for radius in self.n_mols]) # mean number of molecules
        r_renormalized = (n_mean/(number_density_nm3) * 3 / 4 / np.pi * 1E3)**(1/3)  # in A
        self.radii_renormalized = r_renormalized
        #        plt.plot(self.radii, r_renormalized) renormalized r vs normal r
        #        plt.savefig('test.png')
        #target_i is used
        selected_r_renormalized = r_renormalized[target_i]
        # selected_energies are used

        # computing eps
        coef_poly_r = np.polyfit(1.0 / selected_r_renormalized, selected_energies, 1)
        selected_fitted_energies_r = coef_poly_r[0] / selected_r_renormalized + coef_poly_r[1]
        epsilon_r = c / (c + coef_poly_r[0])
        point_a_r = [0, coef_poly_r[1]]
        point_b_r = [10/selected_r_renormalized[-1], selected_energies[-1]]

        self.epsilon_renormalized = epsilon_r
        print("Slope renormalized: ", coef_poly_r[0])
        print("Polarization energy renormalized: ", coef_poly_r[1])
        print("Dielectric permittivity renormalized:", epsilon_r)


        plt.plot(10 / r_renormalized, mean_energies, LineStyle='none', marker='o', color = 'C1', label = 'polarization renormalized', mfc='none')  # whole curve
        plt.plot(10 / selected_r, selected_fitted_energies_r, color='C1')  # ??
        plt.plot([point_a_r[0], point_b_r[0]],[point_a_r[1],point_b_r[1]], label='interpolation renormalized', color='C1', linestyle=':')

        #  auxilary plot
        plt.plot([point_a_25[0], point_b[0]],[point_a_25[1],point_b[1]], label='eps 2.5', linestyle=':', color='grey')
        plt.plot([point_a_30[0], point_b[0]],[point_a_30[1],point_b[1]], label='eps 3.0', linestyle=':', color='grey')
        plt.plot([point_a_35[0], point_b[0]],[point_a_35[1],point_b[1]], label='eps 3.5', linestyle=':', color='grey')

        plt.xlabel('10/R, A')
        plt.ylabel('polarization, eV')
        plt.xlim([10 / b1, 10 / a1])

        ax1 = plt.gca()
        add_inverse_axis(ax1)
        plt.savefig('EPS_FULL_ENV_CUSTOM_NEW.png', dpi=600)
        plt.close()

def plot_4_plus_1(x, Pi, p, y_name, filename, labels_Pi, label_p, xlim):

    for i in range(len(Pi)):
        plt.plot(x, Pi[i], label = labels_Pi[i])
    plt.plot(x, p, label=label_p)
    plt.xlabel('$10/R,  10 / \AA^{-1} $')
    plt.ylabel('{}, eV'.format(y_name))
    plt.legend()
    plt.xlim(xlim)
    plt.grid()
    ax1 = plt.gca()
    add_inverse_axis(ax1)
    plt.savefig(filename)
    plt.close()

def plot_4_plus_1_for_p_vs_n(x, Pi, p, y_name, filename, labels_Pi, label_p, xlim):

    for i in range(len(Pi)):
        plt.plot(x, Pi[i], label = labels_Pi[i])
    plt.plot(x, p, label=label_p)
    plt.xlabel('$N^{-1/3} } $')
    plt.ylabel('{}, eV'.format(y_name))
    plt.legend()
    plt.xlim(xlim)
    plt.grid()
    ax1 = plt.gca()
    add_inverse_axis_n(ax1)
    plt.savefig(filename)
    plt.close()

def plot_1(x, y, filename, y_name, xlim):
    plt.plot(x, y)
    plt.xlabel('$10/R,  10 / \AA^{-1} $')
    plt.ylabel('{}, eV'.format(y_name))
    plt.legend()
    plt.xlim(xlim)
    plt.grid()
    ax1 = plt.gca()
    add_inverse_axis(ax1)
    plt.savefig(filename)
    plt.close()

def add_inverse_axis(initial_axis,
                     rs_plot=np.array([1, 2, 3, 4, 5, 7, 10, 20, 30, 40, 50]),
                     rs_grid=np.array([]),
                     x_label='$R, \AA$'):

    def inv(x):
        return 10 / x

    ax2 = initial_axis.twiny()
    ax2.set_xlabel(x_label)
    ax2.set_xticks(inv(rs_plot), minor=False)
    ax2.set_xticks(inv(rs_grid), minor=True)
    ax2.set_xticklabels(rs_plot)
    ax2.set_xbound(initial_axis.get_xbound())
    ax2.grid(True, which='minor', color = 'blue')
    ax2.grid(True, which='major', color = 'red', linestyle='dotted')

def add_inverse_axis_n(initial_axis, rs_plot=np.array([1, 2, 3, 4, 5, 7, 10, 20, 30, 50, 100]),
                     rs_grid=np.array([])):
    def minus_033(x):
        return x**(-1/3)

    ax2 = initial_axis.twiny()
    ax2.set_xlabel('number of molecules, $N$')
    ax2.set_xticks(minus_033(rs_plot), minor=False)
    ax2.set_xticks(minus_033(rs_grid), minor=True)
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


def load_yaml_from_gzip(file):
    with gzip.open(file, 'rb') as stream:
        try:
            file_content = yaml.load(stream, Loader=yaml.FullLoader)
        except yaml.YAMLError as exc:
            print(exc)
    return file_content

if __name__ == '__main__':
    main()
