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

dis_sim = CurrTempDMRset(rates=rates, analyse_Ecoul=args.Ecoul)
dis_sim.get_nonconv_joblist()