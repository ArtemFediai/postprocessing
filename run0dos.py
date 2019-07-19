from __future__ import absolute_import
import os
import numpy as np
import yaml
import matplotlib
import math

matplotlib.use('Agg')

max_iter = 1000000000
load = False
morph_only = False
collect_current = False
create_plots = False
replay = None
npc = 1500

T = 300
n_cpus = 1
n_r = 30  # 45

num_iter = np.int(1E6)
n = 3  # technical
Ndop = 4+n*3
dop = np.array(np.logspace(np.log10(1E-4), np.log10(1E-1), Ndop))

# a for each doping
Np_I_wish = 200
a = np.array((Np_I_wish / dop) ** (1.0 / 3.0), dtype=np.int)
for i in range(len(a)):
    if a[i] < 20:
        a[i] = 20.0
num_part_after = np.array(a * a * a * dop, dtype=np.float)
np.savetxt('a.txt', a)
np.savetxt('doping.txt', dop)
np.savetxt('n_r_n_dop.txt', [n_r, Ndop])

if __name__ == '__main__':
    n_cpu_vs_dop = np.zeros(Ndop, dtype=np.int)
    N_CPUS = 0
    jobfile = open("joblist", 'w')

    for i_dop in range(0, len(dop)):
        for i_r in range(0, n_r):
            n_cpus = num_part_after[i_dop] // npc + 1
            N_CPUS += n_cpus
            jobfile.write("%i run0dos.run_kmc %i %i\n" % (n_cpus, i_dop, i_r))
            n_cpu_vs_dop[i_dop] = n_cpus
    jobfile.close()
    print('num of doping points = %s' % Ndop)
    print('a = %s' % a)
    print('actual dopant number =  %s' % np.array(num_part_after, dtype = np.int))
    print('doping = %s' % dop)
    print('n_cpus_vs_dop = %s' % n_cpu_vs_dop)
    print('N_CPUS = %s' % N_CPUS)
    print('number of nodes = %s' % (N_CPUS // 20))
    print('num replicas = %s' % n_r)


def run_kmc(i_dop, i_r):
    from lightforge.lightforge import main
    fid = open('no.yml', 'r')
    settings = yaml.load(fid)
    fid.close()
    DirName = 'dop_{}/r_{}'.format(i_dop, i_r)
    if not os.path.isdir(DirName):
        os.makedirs(DirName)
    from shutil import copyfile
    copyfile('analysis_autorun', '{}/analysis_autorun'.format(DirName))

    ###########
    settings["layers"][0]["thickness"] = np.int(a[i_dop])
    settings["morphology_width"] = np.int(a[i_dop])
    settings['layers'][0]['molecule_species'][0]['concentration'] = float(1.0 - dop[i_dop])
    settings['layers'][0]['molecule_species'][1]['concentration'] = float(dop[i_dop]) 
    ###########

    with open('{}/no.yml'.format(DirName), 'w') as fid:
        yaml.dump(settings, fid)
    fid.close()
    os.chdir(DirName)
    main(load=load, replay=replay, settings_file="no.yml", morphology_only=morph_only, num_iterations=num_iter,
         collect_current_density=collect_current, collect_density_plots_only=create_plots)
