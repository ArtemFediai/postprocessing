#!/usr/bin/env python3
import argparse
from SimAnalysis import CurrSysSizeSimulation

# This is an example on how to use the classe CurrSysSizeSimulation from SimAnalysis.
# CurrSysSizeSimulation has to be used in a directory containing folders of
# the format sys_*/r_*. Generates data and plots for each J(l) simulation
# seperately. Shall be used if one wants to analyse the influence of the
# system size on the current for a given temperature, field, DMR and
# disorder.

parser = argparse.ArgumentParser()
parser.add_argument('-marcus', action='store_true', help='Specify if simulation was done with Marcus rates.')
args = parser.parse_args()
print(args)

if args.marcus:
    rates = "Marcus"
else:
    rates = "Miller"
    
simulation = CurrSysSizeSimulation(rates=rates)

simulation.get_av_current()
simulation.plot_av_current(errorbar=True,y_logmin=1,y_logmax=14)
# simulation.plot_av_current(errorbar=True,plot_log=False)

simulation.get_conv_analysis()
simulation.plot_conv_analysis()
