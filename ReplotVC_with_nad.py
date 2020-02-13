"""

Replay of V_C from LF for SID paper

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


    # change here to the folder with data
    folder = 'results/material/energy_levels/DataFigure1/'

    V_C = np.loadtxt(folder + 'V_coul.dat')
    x = np.loadtxt(folder+'dist.dat')
    xxx = np.loadtxt(folder+'xxx.dat')
    nad = np.loadtxt('results/material/energy_levels/id_file_distance_data.dat')[:,3] #  nearest atom distance

    coef = -14.3  # eV/nm

    make_scatter_plot_4VC(
        [[xxx, coef / xxx / 2], [xxx, coef / xxx / 3], [xxx, coef / xxx / 4], [x, V_C], [nad, V_C], [0.5*(x + nad), V_C]],
        "Figure_4_SID_paper.png",
        "scatter",
        labels=["$\epsilon$ = 2", "$\epsilon$ = 3", "$\epsilon$ = 4", "cog", "nad", "0.5x(cog + nad)"],
        xlabel='Pair distance [A]',
        ylabel='V_coulomb [eV]',
        style=['-', '-', '-', '.', '.', 'x'],
        color_dict=['grey', 'C3', 'grey', 'C3', 'C4', 'black'],
        ylim=[-1.5, 0],
        xlim=[0, 30])

# FUNCTIONS

def make_scatter_plot_4VC(data, filename, plottype="graph", errorbars=None, xlabel='Pair distance [nm]',
                          ylabel='|J| [eV]', log=False, labels=np.zeros(30, dtype=str), xlim=None, ylim=None,
                          is_graph=None, color_dict=['C1'] * 8, style=[' '] * 30):  # beliebige anzahl Scatter plots
    # mit  data=[ [x1array,y1array], [x2array,y2array],....]  format also z.b.
    #color_dict=['blue','black', 'C0','C1','C2','C3', 'red', 'C4','C5','C6']
    import matplotlib.pyplot as plt
    if plottype == "graph":
        for i in range(len(data)):
            if errorbars is None:
                plt.plot(data[i][0],data[i][1],style[i],color=color_dict[i],label=labels[i])
            else:
                plt.errorbar(data[i][0],data[i][1],yerr=errorbars[i],color=color_dict[i],label=labels[i])
    else:
        for i in range(len(data)):
            try:
                if is_graph[i]:
                    plt.plot(data[i][0], data[i][1], style[i], color=color_dict[i], label=labels[i])
                else:
                    plt.plot(data[i][0], data[i][1], style[i], color=color_dict[i], label=labels[i])
            except:
                if style[i] == '-':
                    a = 0.5
                else:
                    a = 1.0
                plt.plot(data[i][0], data[i][1], style[i], color=color_dict[i], label=labels[i], alpha=a)
    if xlim is not None:
        plt.xlim(xlim)
    if ylim is not None:
        plt.ylim(ylim)
    if log:
        plt.yscale('log')
    if labels[0] != '':
        plt.legend()
    plt.xlabel(xlabel, fontsize=20)
    plt.xticks(size=15)
    plt.ylabel(ylabel, fontsize=20)
    plt.yticks(size=15)
    plt.tight_layout()
    plt.savefig(filename)
    plt.close()


if __name__ == '__main__':
    main()

