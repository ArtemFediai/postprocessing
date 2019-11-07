#!/usr/bin/env python3
import argparse
from SimAnalysis import CurrTempSimulation

# This is an example on how to use the classe CurrTempSimulation from SimAnalysis.
# CurrTempSimulation has to be used in a directory containing folders of
# the format temp_*/r_*. Generates data and plots for each J(T) simulation
# seperately. Shall be used if one wants to analyse only one J(T) simulation
# (for only one DMR).
# TODO: add parser for Marcus/Miller/E_coul


parser = argparse.ArgumentParser()

parser.add_argument('-marcus', action='store_true', help='Specify if simulation was done with Marcus rates.')
parser.add_argument('-Ecoul', action='store_true')

args = parser.parse_args()

if args.Ecoul:
    rates = "Marcus"  # currently implemented that way in SimAnalysis
else:
    if args.marcus:
        rates = "Marcus" 
    else:
        rates = "Miller"

simulation = CurrTempSimulation(rates=rates,analyse_Ecoul=args.Ecoul, source_dir='',dest_dir='analysis')


simulation.collect_current_data()

simulation.get_av_current()
simulation.plot_av_current(Tlim_low=250, Tlim_high=350,Jlim_low=10000,Jlim_high=40000,errorbar=True)

simulation.get_av_conductivity()
simulation.plot_av_conductivity(errorbar=True)

simulation.get_act_energy(Tlim_low=150,Tlim_high=300)

simulation.get_conv_analysis()
simulation.plot_conv_analysis()