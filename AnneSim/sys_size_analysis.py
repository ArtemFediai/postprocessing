from SimAnalysis import CurrSysSizeSimulation

# This is an example on how to use the classe CurrSysSizeSimulation from SimAnalysis.
# CurrSysSizeSimulation has to be used in a directory containing folders of
# the format sys_*/r_*. Generates data and plots for each J(l) simulation
# seperately. Shall be used if one wants to analyse the influence of the
# system size on the current for a given temperature, field, DMR and
# disorder.

simulation = CurrSysSizeSimulation()

simulation.get_current()
simulation.plot_current(errorbar=True,y_logmin=-15,y_logmax=15)

simulation.get_conv_analysis()
simulation.plot_conv_analysis()
