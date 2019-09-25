from SimAnalysis import CurrFieldDMRset

# This is an example on how to use the class CurrFieldDMRset from SimAnalysis.
# CurrFieldDMRset has to be used in a directory containing folders of
# the format DMR_*/f_*/r_*. Generates data for each J(E) simulation seperately
# but only generates plots for the whole DMR set (all J(E) simulation 
# with different DMRs). Should be used if one is interested in the influence 
# of different DMRs on the J(T) curves etc..


dis_sim = CurrFieldDMRset()
dis_sim.plot_current(errorbar=True)
dis_sim.plot_conv_analysis()
