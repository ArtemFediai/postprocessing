"""
append one yml with second yml
"""
import numpy as np
import yaml
import os
import shutil

def main():


    fn_a = 'settings_ng.yml'
    fn_b = 'append.yml'
    fn_c = fn_a.split('.')[0] + '_old' + '.yml'

    lib_a = get_yml(fn_a)
    lib_b = get_yml(fn_b)

    key = 'homo_lumo_generator'
    value = lib_b[key]

    lib_a['Analysis'][key] = value

    shutil.copy(fn_a,fn_c)

    push_yaml(lib_a, fn_a)

    print('I am done')


def get_yml(qp_settings_file):
    with open(qp_settings_file, 'r') as fid:
        qp_settings = yaml.load(fid, Loader=yaml.FullLoader)
    return qp_settings


def push_yaml(lib_a, qp_settings_file):
    with open(qp_settings_file, 'w') as fid:
        yaml.dump(lib_a, fid, explicit_start=True)

if __name__ == '__main__':
    main()
