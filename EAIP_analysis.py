import numpy as np
import yaml
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
FileName = 'eaip_69b6a5a2aca38ba2f734a6724ea30826_f64d2a34d873ebd530e9c4528c688913_EA_summary'
#FileName = 'eaip_69b6a5a2aca38ba2f734a6724ea30826_f64d2a34d873ebd530e9c4528c688913_IP_summary'
#FileName = 'eaip_f64d2a34d873ebd530e9c4528c688913_69b6a5a2aca38ba2f734a6724ea30826_EA_summary'
#FileName = 'eaip_f64d2a34d873ebd530e9c4528c688913_69b6a5a2aca38ba2f734a6724ea30826_IP_summary'


def main():
    fid = open('{}.yml'.format(FileName), 'r')
    A = yaml.load(fid)
    fid.close()
    n_pairs = len(A['raw_data'])
    dimer_comdist = np.zeros(n_pairs)
    full_env = np.zeros(n_pairs)
    single_delta = np.zeros(n_pairs)

    #
    area = np.pi*3
    colors = (0,0,0)
    #

    for i, str in enumerate(A['raw_data']):
     print(str)
     dimer_comdist[i] = A['raw_data'][str]['dimer_comdist']
     full_env[i] = A['raw_data'][str]['full_env']
     single_delta[i] = A['raw_data'][str]['single_delta']

#1
    plt.figure('Full_env')
    plt.scatter(dimer_comdist, full_env, alpha=0.5)
    plt.title('full_env')
    plt.xlabel('Center-of-the-mass distance [Angstroms]')
    plt.ylabel('Energy "full_env"  [eV]')
    plt.savefig('full_env.png', dpi=600)
    plt.close()

#2
    plt.figure('single_delta')
    plt.scatter(dimer_comdist, single_delta, alpha=0.5)
    plt.title('single_delta')
    plt.xlabel('Center-of-the-mass distance [Angstroms]')
    plt.ylabel('Energy "single_delta"  [eV]')
    plt.savefig('single_delta.png', dpi=600)
    plt.close()

#1+2
    plt.figure('all')
    plt.scatter(dimer_comdist, single_delta, alpha=0.5, label='single_delta', color='C0')
    plt.scatter(dimer_comdist, full_env, alpha=0.5, label='full_env', color='C1')
    plot_trend(dimer_comdist, single_delta, 0)
    plot_trend(dimer_comdist, full_env, 1)
    plt.title('Analysis/EAIP/{}'.format(FileName))
    plt.xlabel('Center-of-the-mass distance [Angstroms]')
    plt.ylabel('Energy "single_delta"  [eV]')
    plt.legend()
    plt.savefig('all_{}.png'.format(FileName), dpi=600)
    plt.close()

def plot_trend(x,y,color_num = 'C0'):
    z = np.polyfit(x, y, 1)
    p = np.poly1d(z)
    plt.plot(x, p(x), color='C{}'.format(color_num), linestyle=':')


if __name__ == '__main__':
    main()