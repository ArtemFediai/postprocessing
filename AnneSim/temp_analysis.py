#!/usr/bin/env python3
import argparse
from SimAnalysis import CurrTempSimulation, CurrTempDMRset

# This is an example on how to use the classe CurrTempSimulation and CurrTempDMRset 
# from SimAnalysis.
#
# CurrTempDMRset has to be used in a directory containing folders of
# the format DMR_*/temp_*/r_*. Generates data for each J(E) simulation seperately
# but only generates plots for the whole DMR set (all J(E) simulation 
# with different DMRs). Should be used if one is interested in the influence 
# of different DMRs on the J(T) curves and E_A vs. DMR.
# Default analysis is done for CurrTemDMRset, Miller rates and no Ecoul consideration.
# Use "-marcus" for Marcus rates. Use "-Ecoul" for Ecoul analysis.
# 
# Use "-single" for CurrTempSimulation.
# CurrTempSimulation has to be used in a directory containing folders of
# the format temp_*/r_*. Generates data and plots for each J(T) simulation
# separately. Shall be used if one wants to analyse only one J(T) simulation
# (for only one DMR). 
# 
# TODO: Caution: make sure to use "-marcus" and "-Ecoul" correctly for "-single", 
#       currently it is not tested if input from "sim_param.txt" matches 
#       the given parsing arguments



parser = argparse.ArgumentParser()

parser.add_argument('-single',action='store_true', help='Specify if one single temperature set is analysed.')
parser.add_argument('-marcus', action='store_true', help='Specify if simulation was done with Marcus rates.')
parser.add_argument('-Ecoul', action='store_true', help='Specify if sim. was done for a specific Ecoul.')


args = parser.parse_args()
print(args)

if args.Ecoul:
    rates = "Marcus"  # currently implemented that way in SimAnalysis
else:
    if args.marcus:
        rates = "Marcus" 
    else:
        rates = "Miller"

if args.single:
    simulation = CurrTempSimulation(rates=rates,analyse_Ecoul=args.Ecoul, source_dir='',dest_dir='analysis')

    simulation.collect_current_data()

    simulation.get_av_current()
    simulation.plot_av_current(Tlim_low=250, Tlim_high=350,Jlim_low=10000,Jlim_high=40000,errorbar=True)

    simulation.get_av_conductivity()
    simulation.plot_av_conductivity(errorbar=True)

    simulation.get_act_energy(Tlim_low=150,Tlim_high=300)

    simulation.get_conv_analysis()
    simulation.plot_conv_analysis()
else:
    try:
        dis_sim = CurrTempDMRset(rates=rates, analyse_Ecoul=args.Ecoul)

        dis_sim.plot_conv_analysis()
        dis_sim.plot_av_current(Tlim_low=250, Tlim_high=400,errorbar=True)

        dis_sim.plot_av_conductivity(errorbar=True,plot_log=True)
        dis_sim.get_act_energy(Tlim_low=250,Tlim_high=350)
        dis_sim.plot_act_energy()
    except:
        print("The format of the simulation does not match. Please check if the right parser arguments for this simulation (-single, -marcus, -Ecoul) were given.")