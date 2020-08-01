from __future__ import print_function
from __future__ import absolute_import
import numpy as np
import matplotlib
from IP_or_EA_vs_R_ver_4 import QPOutput
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
from IP_or_EA_vs_R_ver_4 import Timer
from joblib import Parallel, delayed
import multiprocessing


def main():

    ene_my = np.loadtxt('../aNPD_1_molecule_outershell_more/core_energies.txt')
    ene_my_015 = np.loadtxt('../aNPD_1_mol_015_decay/core_energies.txt')
    ene_my_015_more = np.loadtxt('../aNPD_1_mol_015_decay_more/core_energies.txt')
    plt.plot(ene_my[0] - ene_my[2], label='0.3')
    plt.plot(ene_my_015[0] - ene_my_015[2], label='0.15')
    plt.plot(ene_my_015_more[0] - ene_my_015_more[2], label='0.15 more it.')
    plt.legend()
    plt.xlim([3,30])
    plt.ylim([2, 4.65])
    plt.savefig('test.png')


    with Timer() as t:
        qp_settings_file = 'settings_ng.yml'  # default name of the QP settings file
        with open(qp_settings_file, 'r') as fid:
            qp_settings = yaml.load(fid, Loader=yaml.FullLoader)
        R = int(qp_settings['System']['Shells']['0']['cutoff'])
        ab = [R // 2, R]  # ab is the left/right limit of the epsilon evaluation interval


        Analysis_output = QPOutput()  # ini: returns target folders

    core = Analysis_output.core_ids[0]

    rt_name = 'quantumpatch_runtime_files'
    pmn = ['Mol_'+np.str(core)+'_C_1', 'Mol_'+np.str(core)+'_C_-1', 'uncharged']  #plus minus null
    pmn = [rt_name + '/' + folder for folder in pmn]

    pmn_folders = []

#    energies = np.zeros([3,18]) # array of 3 lists

    len_pmn = len(pmn)

    starter = 0
    for i in range(len_pmn):
        folder = pmn[i]
        pmn[i] = folder
        folders = os.listdir(pmn[i])
        folders_int = np.sort([int(folder) for folder in folders])
        if i == 2:
            folders_int = folders_int[1:]
        if starter == 0:
            energies = np.zeros([3,len(folders)])
            len_folders = len(folders_int)
            starter = 1
        pmn_folders.append(folders_int)
        dummy_e = np.zeros(len_folders)

        print('important:', len_folders)

        ###################################
        # for j in range(len_folders):
        #     subfolder = folders_int[j]
        #     dummy_e[j] = return_energies(folder, subfolder, j)
        ##########################3##########
        ####

        num_cores = multiprocessing.cpu_count()
        j_subfolders = [np.array([j,subfolder]) for [j,subfolder] in np.transpose(np.array([range(len_folders), folders_int]))]


        t0 = time.time()
        j_ene_array = np.array(Parallel(n_jobs=num_cores)(delayed(return_energies)(folder,j_subfolder) for j_subfolder in j_subfolders))
        t1 = time.time() - t0
        print('elapsed time:', (t1-t0))
        ####

        # dummy_e = np.sort(j_ene_array,1)[:,1]
        dummy_e = j_ene_array[:,1]


        energies[i] = np.array(dummy_e)


    np.savetxt('core_energies.txt', energies)
#    plt.show()
    plt.plot(folders_int, energies[0] - energies[2])
    plt.legend()
    plt.xlim([3,20])
    plt.ylim([4.0, 4.8])
    plt.savefig('test_par.png')
#    plt.plot(energies[1][:] - energies[0][:])




    print("\nI am done")
    print("Total Computation Time: {} sec".format(t.interval))


def load_yaml_from_gzip(file):
    with gzip.open(file, 'rb') as stream:
        try:
            file_content = yaml.load(stream, Loader=yaml.FullLoader)
        except yaml.YAMLError as exc:
            print(exc)
    return file_content

def return_energies(folder, j_subfolder):
    j = j_subfolder[0]
    subfolder = j_subfolder[1]
    with Timer() as t1:
        loaded_file = load_yaml_from_gzip(folder + '/' + str(subfolder) + '/' + 'energies.ene.yml.gz')
    print('time:', t1.interval)
    print(subfolder)
    return np.array([j, np.sum(np.array([loaded_file[my_key]["total"] for my_key in loaded_file.keys()]))])




if __name__ == '__main__':
    main()