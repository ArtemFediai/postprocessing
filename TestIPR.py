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

    my_pattern = 'EAIP_for_radius_'

    folders, radii, dict_radii_folder = return_target_folders(my_pattern)

    print("I am done")

    pass

############ FUNCTIONS ################
def return_target_folders(pattern):
    """
    return all folders of the kind "d{2}"+pattern"
    :param: pattern (string)
    :returns: names of folders
    """

    r = re.compile(pattern + '[0-9]*\.[0-9]*')

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



    return match_folders, radii, my_dict


if __name__ == '__main__':
    main()