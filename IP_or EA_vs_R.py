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
matplotlib.use('pdf')
import matplotlib.pyplot as plt
import os
import re
import yaml
import scipy.constants as const
import time
import gzip

def main():


    # my_mol_name = '7b043a1428ff586ba119947b48219969' #TPD
    my_mol_name = 'f273d730b77306e88fc77abb9fa56671' #aNPD
    ab = [10, 20]

    # get ea/ip output
    ip_output = QPOutput(my_mol_name, IP_or_EA='IP')  # ini: returns target folders for IP
    ea_output = QPOutput(my_mol_name, IP_or_EA='EA')  # ini: returns target folders for EA

    # extract averaged energies
    ip_output.extract_mean_energies(mol_name=my_mol_name)
    ea_output.extract_mean_energies(mol_name=my_mol_name)

    # plot ip and ea vs R and 1/R
    ip_output.plot_energies(full_env_or_single_delta='single_delta')
    ip_output.plot_energies(full_env_or_single_delta='full_env')
    ea_output.plot_energies(full_env_or_single_delta='single_delta')
    ea_output.plot_energies(full_env_or_single_delta='full_env')

    # save ip and ea
    ea_output.save()
    ip_output.save()

    # extract and plot eps
    ip_output.extract_eps(ab, full_env_or_single_delta='single_delta')
    ip_output.extract_eps(ab, full_env_or_single_delta='full_env')
    ea_output.extract_eps(ab, full_env_or_single_delta='single_delta')
    ea_output.extract_eps(ab, full_env_or_single_delta='full_env')

    # get vacuum energies from runtime output
    my_runtime_output = RunTimeOutput()
    my_runtime_output.extract_vacuum_energies()  # creates energies_dict with vacuum energies

    # get vacuum IP, EA
    my_runtime_output.compute_vacuum_binding_energies()

    print("\nI am done")


############ FUNCTIONS ################
class QPOutput:
    def __init__(self, mol_name, IP_or_EA='IP'):
        self.mol_name = mol_name
        self.folders = []
        self.IP_or_EA = IP_or_EA
        self.radii = np.array([])
        self.dict_radii_folder = []

        self.return_analysis_folders()

        self.n_mol = len(self.folders)
        self.mean_single_delta = np.zeros(len(self.folders))
        self.mean_full_env = np.zeros(len(self.folders))
        self.binding_energies_dict = \
            {'radii': self.radii, 'single_delta': np.empty(len(self.folders)), 'full_env': np.empty(len(self.folders))}
        self.energies_dict = \
            {'radii': self.radii, 'single_delta': np.empty(len(self.folders)), 'full_env': np.empty(len(self.folders))}


    def return_analysis_folders(self):
        folder_list = os.listdir("Analysis")

        r = re.compile(self.IP_or_EA + '_' + self.mol_name + '_for_radius_' '[0-9]*\.[0-9]*' + '_classical_correction')

        match_folders = [folder for folder in folder_list if r.match(folder)]

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

    def extract_eps(self, ab=[0, 50], full_env_or_single_delta='single_delta'):
        # fit the slope --> get the eps
        if full_env_or_single_delta == "full_env":
            mean_energies = self.mean_full_env  # TODO: all energies to IP, EA
        elif full_env_or_single_delta == "single_delta":
            mean_energies = self.mean_single_delta  # TODO: all energies to IP, EA
        else:
            raise Warning("\"full_env_or_single_delta\" must be either \"full_env\" or \"single_delta\"")

        c = 7.1997176734999995  # coefficient in the formula
        a, b = ab[0], ab[1]

        target_i = np.where(np.logical_and(self.radii > a, self.radii < b))[0]
        selected_r = self.radii[target_i]  # selected interval to compute epsilon
        selected_energies = mean_energies[target_i]  # TODO: all energies to IP/EA

        coef_poly = np.polyfit(1.0 / selected_r, -selected_energies, 1)
        epsilon = (c / (c - coef_poly[0]), c / (c + coef_poly[0])) [full_env_or_single_delta == "single_delta"]
        selected_fitted_energies = coef_poly[0] / selected_r + coef_poly[1]

        print("\nAnalysis of {}(R). Type: {}".format(self.IP_or_EA, full_env_or_single_delta))
        print("Slope: ", coef_poly)
        print("Dielectric permittivity :", epsilon)

        plt.plot(10 / self.radii, -mean_energies, LineStyle='-', marker='o')  # whole curve
        plt.plot(10 / selected_r, selected_fitted_energies)  # ??
        plt.xlabel('10/R, A')
        plt.ylabel('{}, eV'.format(self.IP_or_EA))
        ym_ym = plt.gca().get_ylim()
        plt.plot(10.0/np.array([20, 20]), ym_ym, label='20 A')  # TODO: make it professional
        plt.plot(10.0/np.array([30, 30]), ym_ym, label='30 A')
        plt.plot(10.0/np.array([40, 40]), ym_ym, label='40 A')
        plt.plot(10.0/np.array([50, 50]), ym_ym, label='50 A')
        plt.plot(10.0/np.array([60, 60]), ym_ym, label='60 A')
        plt.legend()
        plt.savefig('{}_vs_1_over_R_{}_epsilon.png'.format(self.IP_or_EA, full_env_or_single_delta))
        plt.close()

    def plot_energies(self, full_env_or_single_delta='single_delta'):
        # fit the slope --> get the eps
        if full_env_or_single_delta == "full_env":
            mean_energies = self.mean_full_env  # TODO: all energies to IP, EA
        elif full_env_or_single_delta == "single_delta":
            mean_energies = self.mean_single_delta  # TODO: all energies to IP, EA
        else:
            raise Warning("\"full_env_or_single_delta\" must be either \"full_env\" or \"single_delta\"")

        x = self.IP_or_EA
        plt.plot(self.radii, -mean_energies, LineStyle='-', marker='o')
        plt.xlabel('R, A')
        plt.ylabel('{}, eV'.format(x))
        plt.savefig('{}_vs_R_{}.png'.format(x, full_env_or_single_delta))
        plt.close()

        plt.plot(10/self.radii, -mean_energies, LineStyle='-', marker='o')
        plt.xlabel('10/R, A')
        plt.ylabel('{}, eV'.format(x))
        plt.savefig('{}_vs_1_over_R_{}.png'.format(x, full_env_or_single_delta))
        plt.close()

    def plot_full_env(self):
        x = self.IP_or_EA
        plt.plot(self.radii, -self.mean_full_env, LineStyle='-', marker='o')
        plt.xlabel('R, A')
        plt.ylabel('{}, eV'.format(x))
        plt.savefig('{}_vs_R_fe.png'.format(x))
        plt.close()

        plt.plot(10/self.radii, -self.mean_full_env, LineStyle='-', marker='o')
        plt.xlabel('10/R, A')
        plt.ylabel('{}, eV'.format(x))
        plt.savefig('{}_vs_1_over_R_fe.png'.format(x))
        plt.close()

    def save(self):
        x = self.IP_or_EA
        print("I save data to {}.dat and R.dat".format(x))
        np.savetxt('r.dat', self.radii)
        np.savetxt('{}_single_delta.dat'.format(x), -self.mean_single_delta)
        np.savetxt('{}_full_env.dat'.format(x), -self.mean_full_env)


class RunTimeOutput:
    def __init__(self):
        self.return_runtime_folders()
        self.vacuum_energies_dict = {'positive': [], 'negative': [], 'neutral': []}  # ini

    def return_runtime_folders(self):
        """
        returns paths to the folders of the last charged/neutral steps
        :return:
        positive folders
        negative folders
        neutral folder
        """
        def append_last_folder(path):
            path = np.atleast_1d(path)
            path_new = []
            for i in range(len(path)):
                step_list = os.listdir(path[i])
                all_numbers = [int(number) for number in step_list]
                number = max(all_numbers)
                item = path[i] + '/' + str(number)
                path_new.append(item)
            return path_new

        path_prefix = 'quantumpatch_runtime_files'
        all_folders = os.listdir(path_prefix)
        pos_paths = [path_prefix + '/' + folder for folder in all_folders if folder.endswith("_C_1")]
        neg_paths = [path_prefix + '/' + folder for folder in all_folders if folder.endswith("_C_-1")]

        pos_paths = append_last_folder(pos_paths)
        neg_paths = append_last_folder(neg_paths)

        if len(pos_paths) == len(neg_paths):
            Warning("number of negative molecules is not the same as the number of positive. Is this what you expected?")

        self.n_mol = len(pos_paths)
        neutral_path = path_prefix + '/uncharged'
        neutral_path = append_last_folder(neutral_path)[0]

        self.positive_folders = pos_paths
        self.negative_folders = neg_paths
        self.neutral_folder = neutral_path

        # print('neutral folder: \n', self.neutral_folder)
        # print('positive folders: \n', self.positive_folders)
        # print('negative folders: \n', self.negative_folders)
    def extract_vacuum_energies(self):

        def extract_charged_energies(folders, yaml_name):
            nested_dict = {}
            for i, folder in enumerate(folders):
                path2yaml = folder + '/' + yaml_name
                energies = load_yaml_from_gzip(path2yaml)
                mol_number = path2yaml.split('/')[1].split('_')[1]
                energy = energies[mol_number]['total']
                nested_dict[mol_number] = energy
            return nested_dict

        def extract_neutral_energies(folder, yaml_name, positive_folders):
            nested_dict = {}
            for i in range(self.n_mol):
                path2yaml = folder + '/' + yaml_name
                energies = load_yaml_from_gzip(path2yaml)
                mol_number = positive_folders[i].split('/')[1].split('_')[1]
                energy = energies[mol_number]['total']
                nested_dict[mol_number] = energy
            return nested_dict

        yaml_name = 'energies.ene.yml.gz'  # yaml file (first-time)

        vacuum_energies_lib_filename = "quantumpatch_runtime_files/vacuum_energies_lib.yaml"
        if not os.path.exists(vacuum_energies_lib_filename):
            print('\nExtracting vacuum energies. Please wait ...')
            # negative
            self.vacuum_energies_dict['negative'] = extract_charged_energies(self.negative_folders, yaml_name)
            print("anions energies extracted")
            # positive
            self.vacuum_energies_dict['positive'] = extract_charged_energies(self.positive_folders, yaml_name)
            print("cations folders extracted")
            # neutral
            self.vacuum_energies_dict['neutral'] = extract_neutral_energies(self.neutral_folder, yaml_name, self.positive_folders)
            print("neutral energies extracted")
            with open(vacuum_energies_lib_filename, 'w+') as fid:
                yaml.dump(self.vacuum_energies_dict, fid, default_flow_style=False)

        else:
            with open(vacuum_energies_lib_filename) as fid:
                self.vacuum_energies_dict = yaml.load(fid, Loader=yaml.FullLoader)
            print("\nQuick load of already extracted vacuum energies")

    def compute_vacuum_binding_energies(self):
        def mean_for_dict(x):
            return np.array(list(x.values())).mean()

        for charged_state in self.vacuum_energies_dict.keys():
            self.vacuum_energies_dict[charged_state]['mean'] = mean_for_dict(self.vacuum_energies_dict[charged_state])

        print('nix')



def load_yaml_from_gzip(file):
    with gzip.open(file, 'rb') as stream:
        try:
            file_content = yaml.load(stream, Loader=yaml.FullLoader)
        except yaml.YAMLError as exc:
            print(exc)
    return file_content



if __name__ == '__main__':
    main()