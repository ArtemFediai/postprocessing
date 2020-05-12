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
        Analysis_output = QPOutput(my_mol_name)  # ini: returns target folders for IP

        # get E0, E+, E- (s-d)
        # get E0, E+, E- (s-d)
        # from Analysis folder

        # extract averaged energies
        Analysis_output.extract_mean_energies(mol_name=my_mol_name)


        # get vacuum energies from runtime output
        my_runtime_output = RunTimeOutput()
        my_runtime_output.extract_vacuum_energies(R)  # creates dicts with vacuum/damped energies

        # get vacuum IP, EA
        my_runtime_output.compute_vacuum_binding_energies()

        # compute and plot ea, ip from runtime
        my_runtime_output.compute_vacuum_ea_ip()

        # compute full env damp
        my_runtime_output.compute_damped_full_env_energies()

        # extract eps full env
        my_runtime_output.extract_eps_from_polarization(ab)

        # combine outputs into a single object
        my_output = Output(my_runtime_output, ip_output=ip_output, ea_output=ea_output)

        # compute and plot polarizations
        my_output.compute_polarization()
        my_output.plot_polarization()
        my_output.save_polarization()

        # pcs output
        pcs_pcs_polarization_path = 'polarization_SD_core_pcs/output.yaml'
        if os.path.exists(pcs_pcs_polarization_path):
            pcs_pcs_output = PcsPcsOutput()
            pcs_pcs_output.compute_polarization()
            pcs_pcs_output.plot_polarization()

    print("\nI am done")
    print("Total Computation Time: {} sec".format(t.interval))


############ FUNCTIONS ################
class QPOutput:
    def __init__(self, mol_name):
        self.mol_name = mol_name
        self.folders_IP,\
        self.radii,\
        self.dict_radii_folder_IP = self.return_analysis_folders(IP_or_EA='IP')
        self.folders_EA, \
        self.radii,\
        self.dict_radii_folder_EA = self.return_analysis_folders(IP_or_EA='EA')
        self.extract_core_ids()
        self.n_mol = len(self.core_ids)

        self.E_sd_plus = {}  # E0 + V(core-env)
        self.E_sd_minus = {}
        self.E_sd_0 = {}

        self.E_sum_env_plus = {}  # E(env) + 2*V(env-env) + V(core-env)
        self.E_sum_env_minus = {}
        self.E_sum_env_0 = {}

        # 4 components of energy -->
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

        self.extract_energies(IP_or_EA='IP')
        self.extract_energies(IP_or_EA='EA')

        # old
        self.mean_single_delta = np.zeros(len(self.folders))
        self.mean_full_env = np.zeros(len(self.folders))
        self.binding_energies_dict = \
            {'radii': self.radii, 'single_delta': np.empty(len(self.folders)), 'full_env': np.empty(len(self.folders))}
        self.energies_dict = \
            {'radii': self.radii, 'single_delta': np.empty(len(self.folders)), 'full_env': np.empty(len(self.folders))}

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

    def extract_energies(self, IP_or_EA):
        """
        neutral, charged energies: sd,fe
        :param IP_or_EA:
        :return:
        """
        # name files to save to
        charge_str = '_plus' if IP_or_EA == 'IP' else '_minus'
        charge_int = '1' if IP_or_EA == 'IP' else '-1'
        E_charged = "Analysis/{}_{}.yaml".format(self.mol_name, charge_str)
        E_0 =  "Analysis/{}_{}.yaml".format(self.mol_name, '_0')

        for i, radius in enumerate(self.radii):
            if IP_or_EA == 'IP':
                self.E_sd_plus[radius] = {}
                self.E_sum_env_plus[radius] = {}
                self.V0_plus[radius] = {}
                self.V0_0[radius] = {}
                self.E_env_plus[radius] = {}
                self.E_env_0[radius] = {}
                self.V_env_plus[radius] = {}
                self.V_env_0[radius] = {}
            else:
                self.E_sd_minus[radius] = {}
                self.E_sum_env_minus[radius] = {}
            self.E_sd_0[radius] = {}
            self.E_sum_env_0[radius] = {}
            for j, core_id in enumerate(self.core_ids):
                if IP_or_EA == 'IP':
                    path2file = \
                        'Analysis/' + self.dict_radii_folder_IP[radius] + '/Matrix-analysis-' + IP_or_EA + '_' \
                        + self.mol_name + '-Mol_' + str(core_id) + '_C_' + charge_int + '.yml'
                else:
                    path2file == \
                        'Analysis/' + self.dict_radii_folder_EA[radius] + '/Matrix-analysis-' + IP_or_EA + '_' \
                        + self.mol_name + '-Mol_' + str(core_id) + '_C_' + charge_int[IP_or_EA] + '.yml'

                with open(path2file) as fid:
                    this_dict = yaml.load(fid, Loader=yaml.SafeLoader)
                if IP_or_EA == 'IP':
                    #sd extraction
                    self.E_sd_plus[radius][core_id] = this_dict[int(core_id)]['total_e_charged']
                    self.E_sd_0[radius][core_id] = this_dict[core_id]['total_e_uncharged']
                    self.E0_plus[core_id] = this_dict[int(core_id)]['total_e_charged_vacuum']
                    self.E0_0[core_id] = this_dict[int(core_id)]['total_e_uncharged_vacuum']

                    #sd processing
                    self.V0_plus[radius][core_id] = self.E_sd_plus[radius][core_id] - self.E0_plus[core_id]
                    self.V0_0[radius][core_id] = self.E_sd_0[radius][core_id] - self.E0_0[core_id]

                    #env extraction
                    env_sub_dict = this_dict[int(core_id)]['environment_molecules']
                    all_env_e = [env_sub_dict[env_id]['total_e_charged'] for env_id in env_sub_dict]
                    self.E_sum_env_plus[radius][int(core_id)] = np.sum(all_env_e) if not all_env_e == [] else 0.0
                    all_env_e = [env_sub_dict[env_id]['total_e_uncharged'] for env_id in env_sub_dict]
                    self.E_sum_env_0[radius][int(core_id)] = np.sum(all_env_e) if not all_env_e == [] else 0.0
                    all_env_e = [env_sub_dict[env_id]['total_e_charged_vacuum'] for env_id in env_sub_dict]
                    self.E_env_plus[radius][int(core_id)] = np.sum(all_env_e) if not all_env_e == [] else 0.0
                    all_env_e = [env_sub_dict[env_id]['total_e_uncharged_vacuum'] for env_id in env_sub_dict]
                    self.E_env_0[radius][int(core_id)] = np.sum(all_env_e) if not all_env_e == [] else 0.0

                    #fe and sd processing
                    self.V_env_plus[radius][core_id] = 0.5 * (self.E_sum_env_plus[radius][core_id]
                                                              - self.E_env_plus[radius][core_id]
                                                              - self.V0_plus[radius][core_id])

                    self.V_env_0[radius][core_id] = 0.5 * (self.E_sum_env_0[radius][core_id]
                                                              - self.E_env_0[radius][core_id]
                                                              - self.V0_0[radius][core_id])

                else:
                    self.E_sd_minus[radius][core_id] = this_dict[int(core_id)]['total_e_charged']
                    self.E_sd_0[radius][core_id] = this_dict[int(core_id)]['total_e_uncharged']
                    self.E0_minus[core_id] = this_dict[int(core_id)]['total_e_charged_vacuum']
                    self.E0_0[core_id] = this_dict[int(core_id)]['total_e_uncharged_vacuum']
                    env_sub_dict = this_dict[int(core_id)]['environment_molecules']
                    all_env_e = [env_sub_dict[env_id]['total_e_charged'] for env_id in env_sub_dict]
                    self.E_sum_env_minus[radius][int(core_id)] = np.sum(all_env_e) if not all_env_e == [] else 0.0
                    all_env_e = [env_sub_dict[env_id]['total_e_uncharged'] for env_id in env_sub_dict]
                    self.E_sum_env_0[radius][int(core_id)] = np.sum(all_env_e) if not all_env_e == [] else 0.0

        with open('Analysis/energies.pkl', 'wb') as fid:
            pickle.dump([self.E0_plus, self.E0_0, self.E0_minus,
                         self.V0_plus, self.V0_0, self.V0_minus,
                         self.V_env_plus, self.V_env_0, self.V_env_minus,
                         self.E_env_plus, self.E_env_0, self.E_env_minus],
                         fid)
        # TODO delete to check ifexist above and do load without doing
        with open('Analysis/energies.pkl', 'rb') as fid:
            [self.E0_plus, self.E0_0, self.E0_minus,
             self.V0_plus, self.V0_0, self.V0_minus,
             self.V_env_plus, self.V_env_0, self.V_env_minus,
             self.E_env_plus, self.E_env_0, self.E_env_minus] \
                = pickle.load(fid)
        print('here now')


    def extract_mean_energies(self, mol_name):
        binding_energies_file = "Analysis/{}_{}.yaml".format(self.IP_or_EA, self.mol_name)
        energies_file = "Analysis/minus_{}_{}.yaml".format(self.IP_or_EA, self.mol_name)

        files_exist = os.path.exists(binding_energies_file) | os.path.exists(energies_file)

        if files_exist:
            with open(binding_energies_file, 'r') as fid:
                self.binding_energies_dict = yaml.load(fid, Loader=yaml.FullLoader)
            with open(energies_file, 'r') as fid:
                self.energies_dict = yaml.load(fid, Loader=yaml.FullLoader)
            self.mean_single_delta = np.array(self.energies_dict['single_delta'])
            self.mean_full_env = np.array(self.energies_dict['full_env'])


        else:
            for i, radius in enumerate(self.radii):
                path2file = \
                    'Analysis/' + self.dict_radii_folder[radius] + '/' + self.IP_or_EA + '_' + mol_name + '_summary.yml'
                with open(path2file) as fid:
                    this_dict = yaml.load(fid, Loader=yaml.SafeLoader)
                self.mean_single_delta[i] = this_dict['mean_single_delta']
                self.mean_full_env[i] = this_dict['mean_full_env']

            # dictionary for ip and ea (binding energies)
            self.binding_energies_dict['radii'] = self.radii.tolist()
            self.binding_energies_dict['single_delta'] = (-self.mean_single_delta).tolist()
            self.binding_energies_dict['full_env'] = (-self.mean_full_env).tolist()
            with open(binding_energies_file, 'w+') as fid:
                yaml.dump(self.binding_energies_dict, fid)

            # dictionary for -ip and -ea (energies)
            self.energies_dict['radii'] = self.radii.tolist()
            self.energies_dict['single_delta'] = (self.mean_single_delta).tolist()
            self.energies_dict['full_env'] = (self.mean_full_env).tolist()
            with open(energies_file, 'w+') as fid:
                yaml.dump(self.energies_dict, fid)

    def save(self):
        x = self.IP_or_EA
        print("I save data to {}.dat and R.dat".format(x))
        np.savetxt('r.dat', self.radii)
        np.savetxt('{}_single_delta.dat'.format(x), -self.mean_single_delta)
        np.savetxt('{}_full_env.dat'.format(x), -self.mean_full_env)


class RunTimeOutput:
    def __init__(self):
        self.positive_folders, \
        self.negative_folders, \
        self.neutral_folder = self.return_runtime_folders(path_prefix='quantumpatch_runtime_files', min_or_max='max')
        self.n_mol = len(self.negative_folders)
        self.positive_folders_vacuum, \
        self.negative_folders_vacuum, \
        self.neutral_folder_vacuum = \
        self.return_runtime_folders(path_prefix='vacuum/quantumpatch_runtime_files', min_or_max='min')  # 'vacuum'
        self.vacuum_energies_dict = {'positive': {}, 'negative': {}, 'neutral': {}}  # ini
        self.damped_energies_dict = {'positive': {}, 'negative': {}, 'neutral': {}, 'radii': {}}  # ini
        self.full_env_dict = {}

    def return_runtime_folders(self, path_prefix, min_or_max):
        """
        returns paths to the folders of the last charged/neutral steps
        :return:
        positive folders
        negative folders
        neutral folder
        (all with the last iteration)
        """
        all_folders = os.listdir(path_prefix)
        pos_paths = [path_prefix + '/' + folder for folder in all_folders if folder.endswith("_C_1")]
        neg_paths = [path_prefix + '/' + folder for folder in all_folders if folder.endswith("_C_-1")]
        neutral_path = path_prefix + '/uncharged'
        pos_paths = append_folder(pos_paths, min_or_max)
        neg_paths = append_folder(neg_paths, min_or_max)
        neutral_path = append_folder(neutral_path, min_or_max)[0]
        if not len(pos_paths) == len(neg_paths):
            warnings.warn("number of +/- molecules are not equal!")
        return pos_paths, neg_paths, neutral_path

    def return_all_TM_molecules(self, R, file_name='structurePBC.cml'):  # TODO extract R as the radius
        """
        uses QP parser
        :return:
        all TM molecules idxs
        all core_mol-env_mols idxs
        """

        from Shredder.Parsers.ParserCommon import parse_system
        self.system = parse_system(file_name)
        num_mols_all = len(self.system.cogs)
        mol_idxs_all = np.linspace(0, num_mols_all - 1, num_mols_all, dtype=int).tolist()
        for folder in self.negative_folders:
            core_id = int(folder.split('/')[1].split('_')[1])  # core molecule number
            dists_core_env = np.linalg.norm(self.system.cogs[core_id] - self.system.cogs,
                                            axis=1).tolist()  # R_core - R_env[i]
            self.full_env_dict[core_id] = \
                {mol_id: {'dist': dists_core_env[mol_id]} for mol_id in mol_idxs_all if dists_core_env[mol_id] <= R}
            # includes core molecule!

    def extract_vacuum_energies(self, R):
        def extract_charged_energies(folders, yaml_name, pos_or_neg):
            for i, folder in enumerate(folders):
                path2yaml = folder + '/' + yaml_name
                energies = load_yaml_from_gzip(path2yaml)
                core_id = path2yaml.split('/')[1].split('_')[1]  # core molecule number
                vacuum_energy = energies[core_id]['damped_energies']['vacuum']  # core vacuum energy
                for mol_id in self.full_env_dict[int(core_id)]:  # full env
                    self.full_env_dict[int(core_id)][mol_id][pos_or_neg] = energies[str(mol_id)]['damped_energies'][
                        'vacuum']
                self.vacuum_energies_dict[pos_or_neg][int(core_id)] = vacuum_energy
                del energies[core_id]['damped_energies']['vacuum']  # remove vacuum energy to avoid confusion
                damped_energies_lib = energies[core_id]['damped_energies']
                self.damped_energies_dict[pos_or_neg][int(core_id)] = [damped_energies_lib[key] for key in
                                                                       damped_energies_lib if
                                                                       key.split('_')[1] == core_id]


        def extract_neutral_energies(folder, yaml_name, positive_folders):
            for i in range(self.n_mol):
                path2yaml = folder + '/' + yaml_name
                energies = load_yaml_from_gzip(path2yaml)
                core_id = positive_folders[i].split('/')[1].split('_')[1]
                energy = energies[core_id]['damped_energies']['vacuum']
                for mol_id in self.full_env_dict[int(core_id)]:  # full env
                    self.full_env_dict[int(core_id)][mol_id]['neutral'] = energies[str(mol_id)]['damped_energies'][
                        'vacuum']
                self.vacuum_energies_dict['neutral'][int(core_id)] = energy  # core vacuum energy
                del energies[core_id]['damped_energies']['vacuum']  # remove vacuum energy to avoid confusion
                damped_energies_lib = energies[core_id]['damped_energies']
                self.damped_energies_dict['neutral'][int(core_id)] = \
                    [damped_energies_lib[key] for key in damped_energies_lib if key.split('_')[1] == core_id]
                self.damped_energies_dict['radii'] = \
                    [float(key.split('_')[3]) for key in damped_energies_lib.keys() if key.split('_')[1] == core_id]

        yaml_name = 'energies.ene.yml.gz'  # yaml file (first-time)

        vacuum_energies_lib_filename = "quantumpatch_runtime_files/vacuum_energies_lib.yaml"
        damped_energies_lib_filename = "quantumpatch_runtime_files/damped_energies_lib.yaml"
        full_env_energies_lib_filename = "quantumpatch_runtime_files/full_env_energies_lib.yaml"

        if not os.path.exists(vacuum_energies_lib_filename):
            print('\nExtracting damped and vacuum energies. Please wait ...')
            print('(R ={} A)'.format(R))
            # env mol ids
            self.return_all_TM_molecules(R)
            print("env mols ids and distances extracted")
            # negative
            extract_charged_energies(self.negative_folders, yaml_name, pos_or_neg='negative')
            print("anions energies extracted")
            # positive
            extract_charged_energies(self.positive_folders, yaml_name, pos_or_neg='positive')
            print("cations folders extracted")
            # neutral
            extract_neutral_energies(self.neutral_folder, yaml_name, self.positive_folders)
            print("neutral energies extracted")

            # mean for vacuum energies
            for charged_state in self.vacuum_energies_dict.keys():
                self.vacuum_energies_dict[charged_state]['mean'] = \
                    mean_for_dict(self.vacuum_energies_dict[charged_state])
            # save vacuum energies
            with open(vacuum_energies_lib_filename, 'w+') as fid:
                yaml.dump(self.vacuum_energies_dict, fid, default_flow_style=False)
            # mean for damped energies
            for charged_state in self.vacuum_energies_dict.keys():
                self.damped_energies_dict[charged_state]['mean'] = \
                    mean_for_dict(self.damped_energies_dict[charged_state])
            # save damped energies
            with open(damped_energies_lib_filename, 'w+') as fid:
                yaml.dump(self.damped_energies_dict, fid, default_flow_style=False)
            # save full_env energies
            with open(full_env_energies_lib_filename, 'w+') as fid:
                yaml.dump(self.full_env_dict, fid, default_flow_style=False)

        else:
            with open(vacuum_energies_lib_filename) as fid:
                self.vacuum_energies_dict = yaml.load(fid, Loader=yaml.FullLoader)
            print("\nQuick load of already extracted vacuum energies")

            with open(damped_energies_lib_filename) as fid:
                self.damped_energies_dict = yaml.load(fid, Loader=yaml.FullLoader)
            print("\nQuick load of already extracted damped energies")

            with open(full_env_energies_lib_filename) as fid:
                self.full_env_dict = yaml.load(fid, Loader=yaml.FullLoader)
            print("\nQuick load of already extracted full env vacuum energies")

    def compute_vacuum_binding_energies(self):  # ip or ea
        self.vacuum_ip = self.vacuum_energies_dict["positive"]['mean'] - self.vacuum_energies_dict["neutral"]['mean']
        self.vacuum_ea = self.vacuum_energies_dict["neutral"]['mean'] - self.vacuum_energies_dict["negative"]['mean']

        print('\nvacuum ip', self.vacuum_ip)
        print('vacuum ea', self.vacuum_ea)

    def compute_vacuum_ea_ip(self):  # short hand function
        """
        polarization from runtime
        :return:
        """
        self.e_plus = np.array(self.damped_energies_dict['positive']['mean'])  # damped -->
        self.e_minus = np.array(self.damped_energies_dict['negative']['mean'])
        self.e_0 = np.array(self.damped_energies_dict['neutral']['mean'])
        self.e_plus_vacuum = self.vacuum_energies_dict["positive"]['mean']  # vacuum -->
        self.e_minus_vacuum = self.vacuum_energies_dict["negative"]['mean']
        self.e_0_vacuum = self.vacuum_energies_dict["neutral"]['mean']

        self.r = np.array(self.damped_energies_dict['radii'])

        self.ip = self.e_plus - self.e_0
        self.ea = self.e_0 - self.e_minus

        self.V = self.e_0 - self.e_0_vacuum

        plt.plot(10 / self.r, self.e_0 - self.e_0_vacuum, label='e_0')
        plt.plot(10 / self.r, self.e_plus - self.e_plus_vacuum, label='eplus')
        plt.plot(10 / self.r, self.e_minus - self.e_minus_vacuum, label='eminus')
        # plt.plot(10/self.r, self.e_plus - self.e_plus_vacuum - self.V, label = 'eplus corrected', marker = 'o', color = 'C1')
        # plt.plot(10/self.r, self.e_minus - self.e_minus_vacuum + self.V, label = 'eminus corrected', marker = 'o', color = 'C2')
        plt.plot(10 / self.r, 0.5 * (self.e_plus - self.e_plus_vacuum + self.e_minus - self.e_minus_vacuum),
                 LineStyle='-.', label='mean')
        plt.xlim([0, 0.5])
        plt.ylim([-1.4, -1])

        plt.legend()
        plt.savefig('energies_runtime.png')
        plt.close()

        self.v_plus = self.e_plus - self.e_plus_vacuum
        self.v_0 = self.e_0 - self.e_0_vacuum
        self.v_minus = self.e_minus - self.e_minus_vacuum

        self.p_cation = -(self.v_plus - self.v_0)
        self.p_anion = self.v_0 - self.v_minus

        plt.plot(10 / self.r, self.p_cation, label='cation')
        plt.plot(10 / self.r, self.p_anion, label='anion')

        # plt.plot(10/self.r, self.p_cation + self.v_0, label='cation corrected', marker='.', color='C0', LineStyle='')
        # plt.plot(10/self.r, self.p_anion - self.v_0, label='anion corrected', marker='.', color='C1', LineStyle='')

        plt.plot(10 / self.r, 0.5 * (self.p_cation + self.p_anion), LineStyle='-.', label='mean')

        plt.legend()

        plt.xlim([0, 2.0])
        plt.ylim([-0.1, 1.4])
        plt.grid()
        plt.savefig('polarization.png')
        plt.close()

        plt.plot(10 / self.r, self.ip)
        plt.savefig('ip_from_runtime.png')
        plt.close()

        plt.plot(10 / self.r, self.ea)
        plt.savefig('ea_from_runtime.png')
        plt.close()

    def compute_damped_full_env_energies(self):
        radii = self.damped_energies_dict['radii']  # I use damped radii by default
        n_r = len(radii)
        n_c = len(self.full_env_dict)
        # vacuum things -->
        ip0 = self.vacuum_energies_dict['positive']['mean'] - self.vacuum_energies_dict['neutral']['mean']
        ea0 = self.vacuum_energies_dict['neutral']['mean'] - self.vacuum_energies_dict['negative']['mean']
        # <--vacuum things
        e_plus_full_env, e_minus_full_env, e_0_full_env, ip_fe, ea_fe, polarization_cation, polarization_anion = \
            np.zeros(n_r), np.zeros(n_r), np.zeros(n_r), np.zeros(n_r), np.zeros(n_r), np.zeros(n_r), np.zeros(n_r)
        for i, r in enumerate(radii):
            for i_core, core_id in enumerate(self.full_env_dict.keys()):
                env_ids = [env_id for env_id in self.full_env_dict[core_id] if
                           self.full_env_dict[core_id][env_id]['dist'] <= r]
                tmp = [self.full_env_dict[core_id][env_id]['positive'] for env_id in self.full_env_dict[core_id] if
                       self.full_env_dict[core_id][env_id]['dist'] <= r]  # energies
                e_plus_full_env[i] += sum(tmp)
                tmp = [self.full_env_dict[core_id][env_id]['negative'] for env_id in self.full_env_dict[core_id] if
                       self.full_env_dict[core_id][env_id]['dist'] <= r]
                e_minus_full_env[i] += sum(tmp)
                tmp = [self.full_env_dict[core_id][env_id]['neutral'] for env_id in self.full_env_dict[core_id] if
                       self.full_env_dict[core_id][env_id]['dist'] <= r]
                e_0_full_env[i] += sum(tmp)
            e_plus_full_env[i] /= n_c
            e_minus_full_env[i] /= n_c
            e_0_full_env[i] /= n_c
            ip_fe[i] = e_plus_full_env[i] - e_0_full_env[i]
            ea_fe[i] = e_0_full_env[i] - e_minus_full_env[i]

            polarization_cation[i] = ip0 - ip_fe[i]
            polarization_anion[i] = ea_fe[i] - ea0

        self.dp_full_env_cation = polarization_cation  # negative. this means that this reduces the polarization!
        self.dp_full_env_anion = polarization_anion  # also negative

        print('I am here')
        inv_r = 10.0 / np.array(radii)
        plt.plot(inv_r, self.dp_full_env_anion, label='dp due to env (-)', marker='.')
        plt.plot(inv_r, self.dp_full_env_cation, label='dp due tu env (+)', marker='.')
        plt.plot(inv_r, self.p_anion, label='single-delta (-)', marker='x')
        plt.plot(inv_r, self.p_cation, label='single-delta (+)', marker='x')
        plt.plot(inv_r, self.p_anion + self.dp_full_env_anion, label='full-env (-)', marker='o')
        plt.plot(inv_r, self.p_cation + self.dp_full_env_cation, label='full-env (+)', marker='o')
        plt.legend()
        plt.xlabel('10/R, 10/A')
        plt.ylabel('polarization, eV')
        plt.xlim(right=1.0)
        plt.savefig('polarization_full_env.png')
        plt.close()

        #  save full_env_polarization
        self.p_full_env_anion = self.p_anion + self.dp_full_env_anion
        self.p_full_env_cation = self.p_cation + self.dp_full_env_cation

    def extract_eps_from_polarization(self, ab, a1b1=[5, 1E100]):
        # fit the slope --> get the eps
        self.radii = np.array(self.damped_energies_dict['radii'])
        mean_energies = 0.5 * (self.p_full_env_cation + self.p_full_env_anion)  # mean polarization energy

        c = 7.1997176734999995  # coefficient in the formula
        a, b, a1, b1 = ab[0], ab[1], a1b1[0], a1b1[1]  # analyze interval; plot interval

        target_i = np.where(np.logical_and(self.radii > a, self.radii < b))[0]
        selected_r = self.radii[target_i]  # selected interval to compute epsilon
        selected_energies = mean_energies[target_i]  # TODO: all energies to IP/EA

        coef_poly = np.polyfit(1.0 / selected_r, selected_energies, 1)

        selected_fitted_energies = coef_poly[0] / selected_r + coef_poly[1]

        epsilon = c / (c + coef_poly[0])

        print("Slope FULL-ENV CUSTOM: ", coef_poly[0])
        print("Dielectric permittivity FULL-ENV CUSTOM :", epsilon)

        plt.plot(10 / self.radii, mean_energies, LineStyle='-', marker='o')  # whole curve
        plt.plot(10 / selected_r, selected_fitted_energies)  # ??
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
        plt.savefig('EPS_FULL_ENV_CUSTOM.png')
        plt.close()

class Output:
    def __init__(self, run_time_output, ip_output, ea_output):
        self.run_time_output = run_time_output
        self.ip_output = ip_output
        self.ea_output = ea_output
        #  single-delta -->
        self.polarization_anion_vs_R = {"single_delta": []}
        self.polarization_cation_vs_R = {"single_delta": []}
        self.polarization_physical_vs_R = {"single_delta": []}
        self.polarization_nonphysical_vs_R = {"single_delta": []}
        #  full_env -->
        self.polarization_anion_vs_R = {"full_env": []}
        self.polarization_cation_vs_R = {"full_env": []}
        self.polarization_physical_vs_R = {"full_env": []}
        self.polarization_nonphysical_vs_R = {"full_env": []}

        self.radii = np.array(self.ip_output.binding_energies_dict["radii"])

    def compute_polarization(self):
        self.polarization_anion_vs_R["single_delta"] = \
            np.array(self.ea_output.binding_energies_dict["single_delta"]) - self.run_time_output.vacuum_ea
        self.polarization_cation_vs_R["single_delta"] = \
            -np.array(self.ip_output.binding_energies_dict["single_delta"]) + self.run_time_output.vacuum_ip
        self.polarization_physical_vs_R["single_delta"] = \
            0.5 * (self.polarization_anion_vs_R["single_delta"] + self.polarization_cation_vs_R["single_delta"])
        self.polarization_nonphysical_vs_R["single_delta"] = \
            0.5 * (self.polarization_cation_vs_R["single_delta"] - self.polarization_anion_vs_R["single_delta"])

        print("\nPolarization computed")

    def plot_polarization(self, single_delta_or_full_env='single_delta', name_prefix='QP'):
        # short-hand:
        x = single_delta_or_full_env
        p_a = self.polarization_anion_vs_R[single_delta_or_full_env]
        p_c = self.polarization_cation_vs_R[single_delta_or_full_env]
        p_phys = self.polarization_physical_vs_R[single_delta_or_full_env]
        p_non_p = self.polarization_nonphysical_vs_R[single_delta_or_full_env]
        r = self.radii
        inv_r = np.array(10.0 / r)

        figure, ax1 = plt.subplots()

        ax1.plot(inv_r, p_c, marker='o', label='cation')
        ax1.plot(inv_r, p_a, marker='o', label='anion')
        ax1.plot(inv_r, p_phys, marker='o', label='physical')
        ax1.plot(inv_r, p_non_p, marker='o', label='nonphysical')

        ax1.set_xlim(left=0.0, right=10.0)
        ax1.set_ylim(bottom=0.0)
        ax1.set_xlabel('10/R, 10/A')
        ax1.set_ylabel('Polarization Energy {}, eV'.format(x))
        ax1.legend()

        add_inverse_axis(initial_axis=ax1)

        figure.tight_layout()
        figure.savefig('{}_P_vs_1_over_R_{}.png'.format(name_prefix, x), dpi=600)

    def save_polarization(self):
        dict2save = {
            'radii': self.radii.tolist(),
            'polarization_anion_vs_R': self.polarization_anion_vs_R["single_delta"].tolist(),
            'polarization_cation_vs_R': self.polarization_cation_vs_R["single_delta"].tolist(),
            'polarization_physical_vs_R': self.polarization_physical_vs_R["single_delta"].tolist(),
            'polarization_nonphysical_vs_R': self.polarization_nonphysical_vs_R["single_delta"].tolist()
        }

        with open('polarization_dict_QP.yaml', 'w+') as fid:
            yaml.dump(dict2save, stream=fid)


class PcsPcsOutput:
    def __init__(self, single_delta_or_full_env='single_delta'):
        self.single_delta_or_full_env = single_delta_or_full_env
        self.polarization_anion_vs_R = {}
        self.polarization_cation_vs_R = {}
        self.polarization_physical_vs_R = {}
        self.polarization_nonphysical_vs_R = {}
        self.radii = []
        self.qp_output = {'radii': [],
                          'polarization_anion_vs_R': [],
                          'polarization_cation_vs_R': [],
                          'polarization_physical_vs_R': [],
                          'polarization_nonphysical_vs_R': []}

    def compute_polarization(self):

        pcs_pcs_file = 'polarization_SD_core_pcs/output.yaml'
        if os.path.exists(pcs_pcs_file):
            with open(pcs_pcs_file) as fid:
                pcs_out_dict = yaml.load(fid, Loader=yaml.FullLoader)
            print("\nLoading pre-computed pcs-pcs polarization energies")

            self.polarization_anion_vs_R["single_delta"] = np.array(pcs_out_dict["polarization_anion_vs_R"])
            self.polarization_cation_vs_R["single_delta"] = np.array(pcs_out_dict['polarization_cation_vs_R'])
            self.polarization_physical_vs_R["single_delta"] = \
                0.5 * (self.polarization_anion_vs_R["single_delta"] + self.polarization_cation_vs_R["single_delta"])
            self.polarization_nonphysical_vs_R["single_delta"] = \
                0.5 * (self.polarization_cation_vs_R["single_delta"] - self.polarization_anion_vs_R["single_delta"])
            self.radii = np.array(pcs_out_dict["radii"])
        else:
            Warning("You tried to load pcs-pcs polarization, but the file {} is not found.".format(pcs_pcs_file))

    def plot_polarization(self, single_delta_or_full_env='single_delta', name_prefix='Pcs_Pcs'):
        # short-hand:
        x = single_delta_or_full_env
        p_a = self.polarization_anion_vs_R[single_delta_or_full_env]
        p_c = self.polarization_cation_vs_R[single_delta_or_full_env]
        p_phys = self.polarization_physical_vs_R[single_delta_or_full_env]
        p_non_p = self.polarization_nonphysical_vs_R[single_delta_or_full_env]
        r = self.radii
        inv_r = np.array(10.0 / r)

        figure, ax1 = plt.subplots(figsize=[6, 6])
        ax1.plot(inv_r, p_c, label='cation')
        ax1.plot(inv_r, p_a, label='anion')
        ax1.plot(inv_r, p_phys, label='physical')
        ax1.plot(inv_r, p_non_p, label='nonphysical')
        ax1.set_xlim(left=0, right=2.0)
        ax1.set_ylim(bottom=-0.1, top=1.4)
        ax1.set_xlabel('10/R, 10/A')
        ax1.set_ylabel('Polarization Energy, eV')
        #

        add_inverse_axis(initial_axis=ax1)

        #  begin: plot pcs-pcs at the same figure if exist
        if os.path.exists('polarization_dict_QP.yaml'):
            with open('polarization_dict_QP.yaml') as fid:
                out_dict_QP = yaml.load(fid, Loader=yaml.FullLoader)
            self.qp_output["radii"] = np.array(out_dict_QP['radii'])
            self.qp_output["polarization_anion_vs_R"] = \
                np.array(out_dict_QP["polarization_anion_vs_R"])
            self.qp_output["polarization_cation_vs_R"] = \
                np.array(out_dict_QP["polarization_cation_vs_R"])
            self.qp_output["polarization_physical_vs_R"] = \
                np.array(out_dict_QP["polarization_physical_vs_R"])
            self.qp_output["polarization_nonphysical_vs_R"] = \
                np.array(out_dict_QP["polarization_nonphysical_vs_R"])
            qp_r = self.qp_output["radii"]
            qp_inv_r = 10 / qp_r
            qp_a = self.qp_output["polarization_anion_vs_R"]
            qp_c = self.qp_output["polarization_cation_vs_R"]
            qp_p = self.qp_output["polarization_physical_vs_R"]
            qp_np = self.qp_output["polarization_nonphysical_vs_R"]

            marker_style = '.'
            ax1.plot(qp_inv_r, qp_c, marker=marker_style, label='QP: cation', color='C0', LineStyle='')
            ax1.plot(qp_inv_r, qp_a, marker=marker_style, label='QP: anion', color='C1', LineStyle='')
            ax1.plot(qp_inv_r, qp_p, marker=marker_style, label='QP: physical', color='C2', LineStyle='')
            ax1.plot(qp_inv_r, qp_np, marker=marker_style, label='QP: nonphysical', color='C3', LineStyle='')
        # end: plot the same pcs-pcs if exists

        ax1.grid()
        ax1.legend()
        figure.tight_layout()
        figure.savefig('{}_P_vs_1_over_R_{}.png'.format(name_prefix, x), dpi=600)


def load_yaml_from_gzip(file):
    with gzip.open(file, 'rb') as stream:
        try:
            file_content = yaml.load(stream, Loader=yaml.FullLoader)
        except yaml.YAMLError as exc:
            print(exc)
    return file_content


def add_inverse_axis(initial_axis, rs_plot=np.array([1, 2, 3, 4, 5, 7, 10, 20, 50]),
                     rs_grid=np.array([1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 20, 30, 40, 50, 60, 70])):
    def inv(x):
        return 10 / x

    ax2 = initial_axis.twiny()
    ax2.set_xlabel('R, A')
    ax2.set_xticks(inv(rs_plot), minor=False)
    ax2.set_xticks(inv(rs_grid), minor=True)
    ax2.set_xticklabels(rs_plot)
    ax2.set_xbound(initial_axis.get_xbound())
    ax2.grid(True, which='minor')
    ax2.grid(True, which='major')


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


def mean_for_dict(x):
    return np.array(list(x.values())).mean(0).tolist()

#TRASH


def Retriever(bar):
    return (list(globals().keys()))[list(map(lambda x: id(x), list(globals().values()))).index(id(bar))]

def save_vars_by_their_names(prefix, suffix, *args):
    for arg in args:
        retrieved_thing = mod_retrieve_name(arg)
        retrieved_thing1 = retrieve_name(arg)
        retrieved_thing2 = mod_mod_retrieve_name(arg)
        file_name = prefix + mod_retrieve_name(arg) + suffix
        with open(file_name,'w+') as fid:
            yaml.dump(arg, fid)

def mod_mod_retrieve_name(var):
    callers_local_vars = inspect.currentframe().f_back.f_back.f_back.f_locals.items()
    return [var_name for var_name, var_val in callers_local_vars if var_val is var]

def mod_retrieve_name(var):
    callers_local_vars = inspect.currentframe().f_back.f_back.f_locals.items()
    return [var_name for var_name, var_val in callers_local_vars if var_val is var]

def retrieve_name(var):
    callers_local_vars = inspect.currentframe().f_back.f_locals.items()
    return [var_name for var_name, var_val in callers_local_vars if var_val is var]

if __name__ == '__main__':
    main()
