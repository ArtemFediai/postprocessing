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

address = '/home/bh5670/Desktop/fh2_ws/2019/112019/13112019/aNPDF4TCNQ_diluted/dis_j/quantumpatch_runtime_files/uncharged'

step = 8

fid = open(address+"/"+str(step)+"/"+"energies.ene.yml", 'r')
A = yaml.load(fid)
fid.close()

a = []
a = A.keys()
a = list(A)
N = np.size(a)


print("I am done")
