# joint 1 2 3

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

    CT_dir = 'CT'
    IPEA_dir = 'IPEA'
    analysis_dir = 'Analysis'

    dirs = [CT_dir, IPEA_dir]
    numbers = ['1', '2', '3']
    id1 = '69b6a5a2aca38ba2f734a6724ea30826'
    id2 = 'f64d2a34d873ebd530e9c4528c688913'
    mean_full_env = 'mean_full_env'
    std_full_env = 'std_full_env'

    dic_CT_123 = np.array([{}, {}, {}])
    dic_IP_123 = np.array([{}, {}, {}])
    dic_EA_123 = np.array([{}, {}, {}])

    dic_CT_joint = {}
    dic_IP_joint = {}
    dic_EA_joint = {}

    #CT
    CT_file_to_read_template, CT_file_to_write = create_CT_folder(CT_dir, analysis_dir, id1, id2)

    # IP
    IP_file_to_read_template, EA_file_to_read_template, IP_file_to_write, EA_file_to_write =\
        create_IPEA_folders(IPEA_dir,analysis_dir,id1,id2)

    for i, number in enumerate(numbers):
        CT_file_to_read =  CT_file_to_read_template.format(number)
        dic_CT_123[i] = load_yaml(CT_file_to_read)

        IP_file_to_read = IP_file_to_read_template.format(number)
        EA_file_to_read = EA_file_to_read_template.format(number)

        dic_IP_123[i] = load_yaml(IP_file_to_read)
        dic_EA_123[i] = load_yaml(EA_file_to_read)

    # CT begins
    dic_CT_joint[mean_full_env] = np.mean(np.array([dic_CT_123[i][mean_full_env] for i,_ in enumerate(numbers)])).item()
    dic_CT_joint[std_full_env] = np.mean(np.array([dic_CT_123[i][std_full_env] for i,_ in enumerate(numbers)])).item()
    dic_CT_joint['raw_data'] = {}
    raw_data_dict = {}
    for i, number in enumerate(numbers):
        dic_CT_joint['raw_data'].update(dic_CT_123[i]['raw_data'])
    # CT ends

    # IP begins
    dic_IP_joint[mean_full_env] = np.mean(np.array([dic_IP_123[i][mean_full_env] for i,_ in enumerate(numbers)])).item()
    dic_EA_joint[std_full_env] = np.mean(np.array([dic_EA_123[i][std_full_env] for i,_ in enumerate(numbers)])).item()
    dic_IP_joint['raw_data'] = {}
    dic_EA_joint['raw_data'] = {}
    raw_data_dict = {}
    for i, number in enumerate(numbers):
        dic_IP_joint['raw_data'].update(dic_IP_123[i]['raw_data'])
        dic_EA_joint['raw_data'].update(dic_EA_123[i]['raw_data'])

    push_yaml(IP_file_to_write, dic_IP_joint)
    push_yaml(EA_file_to_write, dic_EA_joint)
    # IP ends


########################################################################################################################
def load_yaml(file_name):
    with open(file=file_name) as fid:
        dict = yaml.load(fid, Loader=yaml.FullLoader)
    return dict

def push_yaml(file_name, dict_name):
    with open(file_name, 'w') as fid:
        yaml.dump(dict_name, fid)
    return dict

def create_IPEA_folders(IPEA_dir, analysis_dir, id1,id2):
    IP_filename = 'IP_{}_summary.yml'.format(id2)
    EA_filename = 'EA_{}_summary.yml'.format(id1)

    IP_file_to_read_template = IPEA_dir + '/' + '{}' + '/' + analysis_dir + '/IP/' + IP_filename
    EA_file_to_read_template = IPEA_dir + '/' + '{}' + '/' + analysis_dir + '/EA/' + EA_filename

    joint_dir_IPEA_name = 'IPEA_joint'

    if not os.path.exists(joint_dir_IPEA_name):
        os.mkdir(joint_dir_IPEA_name)
    if not os.path.exists(joint_dir_IPEA_name + '/' + analysis_dir):
        os.mkdir(joint_dir_IPEA_name + '/' + analysis_dir)
    if not os.path.exists(joint_dir_IPEA_name + '/' + analysis_dir + '/EA'):
        os.mkdir(joint_dir_IPEA_name + '/' + analysis_dir + '/EA')
        os.mkdir(joint_dir_IPEA_name + '/' + analysis_dir + '/IP')

    IP_file_to_write = joint_dir_IPEA_name + '/' + analysis_dir + '/IP/' + IP_filename
    EA_file_to_write = joint_dir_IPEA_name + '/' + analysis_dir + '/EA/' + EA_filename

    return IP_file_to_read_template, EA_file_to_read_template, IP_file_to_write, EA_file_to_write

def create_CT_folder(CT_dir, analysis_dir, id1, id2):
    CT_filename = 'CT_{}_-1,{}_1_summary.yml'.format(id1,id2)
    CT_file_to_read_template =  CT_dir + '/' + '{}' + '/' + analysis_dir + '/CT/' + CT_filename
    joint_dir_CT_name = 'CT_joint'

    if not os.path.exists(joint_dir_CT_name):
        os.mkdir(joint_dir_CT_name)
    if not os.path.exists(joint_dir_CT_name + '/' + analysis_dir):
        os.mkdir(joint_dir_CT_name + '/' + analysis_dir)
    if not os.path.exists(joint_dir_CT_name + '/' + analysis_dir + '/CT'):
        os.mkdir(joint_dir_CT_name + '/' + analysis_dir + '/CT')

    CT_file_to_write = joint_dir_CT_name + '/' + analysis_dir + '/CT/' + CT_filename
    return CT_file_to_read_template, CT_file_to_write

if __name__ == '__main__':
    main()