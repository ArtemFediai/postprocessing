#!/usr/bin/env python3
from SimAnalysis import CurrTempSimulation

# This is an example on how to use the classe CurrTempSimulation from SimAnalysis.
# CurrTempSimulation has to be used in a directory containing folders of
# the format temp_*/r_*. Generates data and plots for each J(T) simulation
# seperately. Shall be used if one wants to analyse only one J(T) simulation
# (for only one DMR).



simulation = CurrTempSimulation(source_dir='',dest_dir='analysis/')

simulation.get_current()
simulation.plot_current(Tlim_low=250, Tlim_high=350,Jlim_low=10000,Jlim_high=40000,errorbar=True)

simulation.get_conductivity()
# simulation.plot_conductivity(errorbar=True)

# simulation.get_act_energy(tlim_low=150,tlim_high=300)

simulation.get_conv_analysis()
# simulation.plot_conv_analysis()