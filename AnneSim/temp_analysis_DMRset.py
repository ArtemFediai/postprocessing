#!/usr/bin/env python3
from SimAnalysis import CurrTempDMRset

# This is an example on how to use the class CurrTempDMRset from SimAnalysis.
# CurrTempDMRset has to be used in a directory containing folders of
# the format DMR_*/temp_*/r_*. Generates data for each J(E) simulation seperately
# but only generates plots for the whole DMR set (all J(E) simulation 
# with different DMRs). Should be used if one is interested in the influence 
# of different DMRs on the J(T) curves and E_A vs. DMR.


dis_sim = CurrTempDMRset()
dis_sim.plot_current(errorbar=True)
dis_sim.plot_conductivity(errorbar=True,plot_log=True)
dis_sim.plot_act_energy()
dis_sim.plot_conv_analysis()