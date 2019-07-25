"""
To be used in the folder where you compute something.
Compare many-body and stupid DOS
"""
from __future__ import print_function
from __future__ import absolute_import
import extract
import numpy as np
import matplotlib
matplotlib.use('pdf')
import matplotlib.pyplot as plt
import os
from shutil import copyfile

logscale = False

dop = np.loadtxt("postprocessing/doping.dat")
n_dop = len(dop)
NE = 1000
energy = np.loadtxt("postprocessing/energy.dat")

HOMO = np.zeros([n_dop, NE])

HOMO = np.loadtxt("postprocessing/HOMO.dat")
LUMO = np.loadtxt("postprocessing/LUMO.dat")
stupid_dos = np.loadtxt("postprocessing/stupid_dos.dat")


for i in range(0, n_dop):
    plt.figure(figsize=[4.5, 3])
    plt.plot(energy, HOMO[i, :], label='HOMO (many-body)')
    plt.plot(energy, LUMO[i, :], label='LUMO (many-body)')
    plt.plot(energy, stupid_dos[i, :], label='DOS (mean-field)')
    plt.legend()
    plt.xlim(bottom=0.0)
    if logscale == True:
        plt.yscale('log')
        plt.ylim(bottom=1E-5)
        plt.ylim(top=1E-1)
    plt.xlim([-6.0, -4.5])
    plt.ylabel("DOS [a.u.]")
    plt.xlabel("Energy [eV]")
    plt.tight_layout()
    plt.savefig("st_dos_comp_{}.png".format(i), dpi = 600)
    plt.close()
