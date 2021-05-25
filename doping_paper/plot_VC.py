import os.path
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import sklearn.linear_model
from collections import namedtuple

import yaml

print("this is plot V_C vs distance")

DopedMaterial = namedtuple('DopedMaterial', 'host host_uid dopant dopant_uid eps_r')
NPD = DopedMaterial('NPD', 'f273d730b77306e88fc77abb9fa56671', 'F4TCNQ', 'ac85c387eec36b98aaab7138b62a8b88', '2.76')
SpiroTAD = DopedMaterial('Spiro-TAD', 'bd683f0bdfe2718f3af430061821f909', 'F6TCNQ', '69b6a5a2aca38ba2f734a6724ea30826', '2.64')
SpiroOMeTAD = DopedMaterial('Spiro-OMeTAD', 'f64d2a34d873ebd530e9c4528c688913', 'F6TCNQ', '69b6a5a2aca38ba2f734a6724ea30826', '2.64')

materials = [NPD, SpiroTAD, SpiroOMeTAD]

EPS_0 = 8.85418782E-12
Q0 = 1.60217662E-19

def main():
    fig=plt.figure()
    gs = fig.add_gridspec(2, hspace=0.1)
    (top_ax, bottom_ax) = gs.subplots(sharex=True)
    fig.add_gridspec(hspace=0)

    is_inverted = False

    for i, material in enumerate(reversed(materials)):
        plot_vc_for_material(material, i, bottom_ax, is_inverted=is_inverted)
        plot_rdf_for_material(material, i, top_ax, is_inverted=is_inverted)

    plt.show()
    plt.close(fig)



    fig1, ax1 = plt.subplots()

def plot_rdf_for_material(material, i, ax, is_inverted=False):
    hatches = ['///', '\\\\\\', '|||', '-', '+', 'x', 'o', 'O', '.', '*']
    path_to_yaml = f'/home/artem/Desktop/PaperDopingFigs/Morphology/{material.host}/morphology_analysis/data.yml'
    with open(path_to_yaml, 'r') as rdf_data:
        whole_dict = yaml.load(rdf_data, Loader=yaml.SafeLoader)
        grid = whole_dict['grid']
        grid = [0.5 * (grid[i] + grid[i + 1]) for i in range(len(grid) - 1)]
        rdf = np.array(whole_dict[f'{material.dopant_uid}_{material.host_uid}']['smooth_rdf'])*1.05
        r1 = whole_dict[f'{material.dopant_uid}_{material.host_uid}']['one_mol_distance']

        if not is_inverted:
            ax.plot(grid, rdf, color=f'C{i}')
            #ax.plot([0.0, 30.0], [1.0, 1.0], color='grey', linestyle=':')
            print(r1)
            grid_fill = [x for x in grid if x <= r1]
            grid_not_fill = np.sort(list(set(grid) - set(grid_fill)))
            ax.fill_between(grid_fill,
                            rdf[0:len(grid_fill)],
                            alpha=1,
                            color='w',
                            hatch=hatches[i],
                            edgecolor=f'C{i}',
                            linewidth=2)
            ax.plot(grid_not_fill, rdf[len(grid_fill)-1:-1], color=f'C{i}')
            #rdf_at_r1 = rdf[len(grid_fill)]
            #ax.plot([r1, r1], [0, rdf_at_r1], linestyle=':')
        else:
            print('not implemented')
            ax.plot(1.0/np.array(grid), rdf, color=f'C{i}')
            ax.plot([0.0, 1/2.5], [1.0, 1.0], color='grey', linestyle=':')



def plot_vc_for_material(material, i, ax, is_inverted=False):
    host, host_id = material.host, material.host_uid
    dopant, dopant_id = material.dopant, material.dopant_uid
    path_to_VC = f'data/Coulomb_binding_energy/{host}/V_coulomb_src.dat'
    path_to_VC_new = f'data/Coulomb_binding_energy/{host}/V_updated.dat'
    path_to_expanded_points = f'data/Coulomb_binding_energy/{host}/V_coulomb_exp_host_dopant.dat'
    exp_points = pd.read_csv(filepath_or_buffer=path_to_expanded_points)
    df_VC = pd.read_csv(filepath_or_buffer=path_to_VC_new)

    pair_distance_type = 'pair cog distance kmc [A]'
    #pair_distance_type = ' pair nad distance [A]'

    df_VC = df_VC.sort_values(by=pair_distance_type)

    #
    # df_new_VC = pd.read_csv(filepath_or_buffer=path_to_VC_new)
    # print(df_new_VC)
    # nad = df_new_VC[' pair nad distance [A]']
    # nad
    # print(nad)
    # print('here')
    #

    # print(df_VC)
    if not is_inverted:
        df_VC.plot(ax=ax,
                   x=pair_distance_type,
                   y=' total Vc [eV]',
                   kind='scatter',
                   color=f'C{i}',
                   label=f'{host}:{dopant}',
                   alpha=0.5,
                   s=8)
        print(material.eps_r)
        x_A, y = coulomb_low(eps_r=np.float(material.eps_r))
        ax.plot(x_A, y, color=f'C{i}')
        ax.set_xlim(0, 30)
    else:
        print('not implemented')
        df_VC.plot(ax=ax,
                   x=pair_distance_type,
                   y=' total Vc [eV]',
                   kind='scatter',
                   color=f'C{i}',
                   label=f'{host}:{dopant}',
                   alpha=0.5,
                   s=8)
        print(material.eps_r)
        x_A, y = coulomb_low(eps_r=np.float(material.eps_r))
        ax.plot(x_A, y, color=f'C{i}')
        ax.set_xlim(0, 30)

    ax.set_ylim(-1, 0)


def coulomb_low(xlim_A=(0,30), eps_r=2.7, n_points=1000):
    x_A = np.linspace(*xlim_A, n_points)
    x = x_A * 1.E-10  #  CI
    if xlim_A[0] == 0.0:
        x[0] = x[0] + np.finfo(np.float).eps/max(x)
    y = -Q0/(4*np.pi*EPS_0*eps_r*x)

    return x_A, y

if __name__ == '__main__':
    main()
