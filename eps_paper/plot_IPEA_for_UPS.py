import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
                               AutoMinorLocator)
from matplotlib.lines import Line2D

# def vac_to_binding(x, workfunction):
#     return x+workfunction
#
# def vac_to_binding_inv(x, workfunction):
#     return x-4.82

def plot_theory_data_binding_energy(name, IP_surface, IP_bulk, EA_surface, EA_bulk, IP_bar, EA_bar,IP_vac, EA_vac, binding_lims, workfunction, ylim=[0, 10]):
    fig, ax = plt.subplots(figsize=(10,5))
    ax2 = ax.twiny()
    binding_energy_lims = binding_lims
    vac_energy_lims = np.array(binding_energy_lims) - workfunction

    IP_plt, = ax.bar(IP_bar[0], height=10, width=IP_bar[1], label="IP", color="C0", alpha=0.6)
    EA_plt, = ax.bar(EA_bar[0], 10, width=EA_bar[1], label="EA", color="C1", alpha=0.6)
    ax.vlines([IP_surface,EA_surface], ylim[0], ylim[1], colors="black", linestyles="dashdot")
    ax.vlines([IP_bulk,EA_bulk], ylim[0], ylim[1], colors="black", linestyles="dotted")
    #vac
    ax.vlines([IP_vac, IP_vac], ylim[0], ylim[1], colors="grey", linestyles="solid")
    ax.vlines([EA_vac, EA_vac], ylim[0], ylim[1], colors="grey", linestyles="solid")

    #secax = ax.secondary_xaxis('top', functions=(vac_to_binding, vac_to_binding_inv))
    #secax.minorticks_on()
    ax.xaxis.set_major_locator(MultipleLocator(1))
    ax.xaxis.set_major_formatter(FormatStrFormatter('%d'))
    ax.xaxis.set_minor_locator(MultipleLocator(0.2))
    ax2.xaxis.set_major_locator(MultipleLocator(1))
    ax2.xaxis.set_major_formatter(FormatStrFormatter('%d'))
    ax2.xaxis.set_minor_locator(MultipleLocator(0.2))

    plt.ylim(ylim)
    ax.set_yticks([])
    ax.set_xlim(vac_energy_lims)
    ax2.set_xlim(binding_energy_lims)
#    secax.set_xlabel("Binding energy [eV]")
    ax.set_xlabel("Energy w.r.t. the vacuum level/eV")
    ax2.set_xlabel("Binding energy/eV")
    ax.set_ylabel("Intensity/a.u.")
    custom_lines = [Line2D([0], [0], color="black", lw=2), IP_plt, EA_plt, Line2D([0], [0], linestyle="dashdot", color="black", lw=2),
                    Line2D([0], [0], linestyle="dotted", color="black", lw=2),Line2D([0], [0], linestyle="solid", color="grey", lw=2)]
    if name=="NPD":
        plt.legend(custom_lines, ["UPS and IPES","IP", "EA", "Surface", "Bulk", "Vacuum"], loc="upper left") #
    else:
        plt.legend(custom_lines, ["UPS and IPES", "IP", "EA", "Surface", "Bulk", "Vacuum"])  #
    #plt.legend()
    plt.savefig("{}.png".format(name))
  #  plt.show()

def plot_theory_data_wrt_vacuum(name, IP_surface, IP_bulk, EA_surface, EA_bulk, IP_bar, EA_bar,IP_vac, EA_vac, vac_lims, ylim=[0, 10]):
    fig, ax = plt.subplots(figsize=(10,5), dpi=150)
    vac_energy_lims = vac_lims

    IP_plt, = ax.bar(IP_bar[0], height=10, width=IP_bar[1], label="IP", color="C0", alpha=0.6)
    EA_plt, = ax.bar(EA_bar[0], 10, width=EA_bar[1], label="EA", color="C1", alpha=0.6)
    ax.vlines([IP_surface,EA_surface], ylim[0], ylim[1], colors="black", linestyles="dashdot")
    ax.vlines([IP_bulk,EA_bulk], ylim[0], ylim[1], colors="black", linestyles="dotted")
    #vac
    ax.vlines([IP_vac, IP_vac], ylim[0], ylim[1], colors="grey", linestyles="solid")
    ax.vlines([EA_vac, EA_vac], ylim[0], ylim[1], colors="grey", linestyles="solid")

    #secax = ax.secondary_xaxis('top', functions=(vac_to_binding, vac_to_binding_inv))
    #secax.minorticks_on()
    ax.xaxis.set_major_locator(MultipleLocator(1))
    ax.xaxis.set_major_formatter(FormatStrFormatter('%d'))
    ax.xaxis.set_minor_locator(MultipleLocator(0.2))

    plt.ylim(ylim)
    ax.set_yticks([])
    ax.set_xlim(vac_energy_lims)
#    secax.set_xlabel("Binding energy [eV]")
    ax.set_xlabel("Energy w.r.t. the vacuum level/eV")
    ax.set_ylabel("Intensity/a.u.")
    custom_lines = [Line2D([0], [0], color="black", lw=2), IP_plt, EA_plt, Line2D([0], [0], linestyle="dashdot", color="black", lw=2),
                    Line2D([0], [0], linestyle="dotted", color="black", lw=2),Line2D([0], [0], linestyle="solid", color="grey", lw=2)]
    plt.legend(custom_lines, ["UPS and IPES","IP", "EA", "Surface", "Bulk", "Vacuum"], loc="upper center") #
    #plt.legend()
    plt.savefig("{}.png".format(name))
    #plt.show()


if __name__=="__main__":
  # aNPD
    IP_vac = -6.016
    EA_vac = -0.239
    IP_surface, IP_bulk = -5.855356585736369, -5.482817320152966 # from 5 cores
    EA_surface, EA_bulk = -0.5487299462795663, -0.7850648317545766
    IP_bar = [np.mean([IP_surface, IP_bulk]), np.abs(IP_bulk - IP_surface) + 2*0.08]  # 1 by bulk surface, 2 by std
    EA_bar = [np.mean([EA_surface, EA_bulk]), np.abs(EA_bulk - EA_surface) + 2*0.06]  # same
    plot_theory_data_binding_energy("NPD",IP_surface, IP_bulk, EA_surface, EA_bulk, IP_bar,EA_bar,IP_vac,EA_vac,binding_lims=[-10,8], workfunction=4.82) # workfunction for Silicon 100

    #C60
    IP_vac = -7.649
    EA_vac = -2.512
    IP_surface, IP_bulk = -7.1592219080499575, -6.621296915771398
    EA_surface, EA_bulk = -3.0894178277449016, -3.6207609752891408
    IP_bar= [np.mean([IP_surface,IP_bulk]), np.abs(IP_bulk-IP_surface) + 2*0.01] # 1 by bulk surface, 2 by std
    EA_bar= [np.mean([EA_surface,EA_bulk]), np.abs(EA_bulk-EA_surface) + 2*0.01] # same
    plot_theory_data_binding_energy("C60",IP_surface, IP_bulk, EA_surface, EA_bulk, IP_bar,EA_bar,IP_vac,EA_vac,binding_lims=[-5.3, 4], workfunction=4.294142857)
   # workfunction computed by secondary electron cutoff and photon energy

    # TCTA
    IP_vac = -6.414
    EA_vac = -0.002
    IP_surface, IP_bulk = -6.17, -5.88
    EA_surface, EA_bulk = -0.27, -0.54
    IP_bar = [np.mean([IP_surface, IP_bulk]), np.abs(IP_bulk - IP_surface) + 2 * 0.05]  # 1 by bulk surface, 2 by std
    EA_bar = [np.mean([EA_surface, EA_bulk]), np.abs(EA_bulk - EA_surface) + 2 * 0.02]  # same
    plot_theory_data_wrt_vacuum("TCTA",IP_surface, IP_bulk, EA_surface, EA_bulk, IP_bar, EA_bar,IP_vac,EA_vac, vac_lims=[-9, 3])