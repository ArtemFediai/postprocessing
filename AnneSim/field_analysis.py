#!/usr/bin/env python3
from SimAnalysis import CurrFieldSimulation

# This is an example on how to use the classe CurrFieldSimulation from SimAnalysis.
# CurrFieldSimulation has to be used in a directory containing folders of
# the format f_*/r_*. Generates data and plots for each J(E) simulation
# seperately. Shall be used if one wants to analyse only one J(E) simulation
# (for only one DMR).

simulation = CurrFieldSimulation(source_dir='',dest_dir='analysis/')

simulation.get_current()
simulation.plot_current(errorbar=True)

simulation.get_conv_analysis()
simulation.plot_conv_analysis()