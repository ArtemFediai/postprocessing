from joblib import Parallel, delayed
import multiprocessing
import os
import numpy as np

inputs = range(10)

folders = os.listdir()


def processInput(folder):
#    a = i
#    print(a)
#    b = os.listdir()
#    print(b)
    c = np.loadtxt(folder + '/' + folder + '.txt').item()
    d = {}
    d[folder] = c
    return d




num_cores = multiprocessing.cpu_count()

results = Parallel(n_jobs=num_cores)(delayed(processInput)(folder) for folder in folders)

pass
# for i, folder in enumerate(folders):
#     print('in {} folder I have found {}'.format(folder, results[folder]))