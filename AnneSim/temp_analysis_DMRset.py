#!/usr/bin/env python3
import argparse
from SimAnalysis import CurrTempDMRset

# This is an example on how to use the class CurrTempDMRset from SimAnalysis.
# CurrTempDMRset has to be used in a directory containing folders of
# the format DMR_*/temp_*/r_*. Generates data for each J(E) simulation seperately
# but only generates plots for the whole DMR set (all J(E) simulation 
# with different DMRs). Should be used if one is interested in the influence 
# of different DMRs on the J(T) curves and E_A vs. DMR.

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


dis_sim = CurrTempDMRset(rates=rates, analyse_Ecoul=args.Ecoul)

dis_sim.plot_conv_analysis()
dis_sim.plot_av_current(Tlim_low=250, Tlim_high=400,errorbar=True)

dis_sim.plot_av_conductivity(errorbar=True,plot_log=True)
dis_sim.get_act_energy(Tlim_low=250,Tlim_high=350)
dis_sim.plot_act_energy()