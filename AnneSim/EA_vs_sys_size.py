#!/usr/bin/env python3
import argparse
from SimAnalysis import CurrTempSimulation, CurrTempDMRset
import glob



parser = argparse.ArgumentParser()
parser.add_argument('-marcus', action='store_true', help='Specify if simulation was done with Marcus rates.')
args = parser.parse_args()
print(args)

if args.marcus:
    rates = "Marcus" 
else:
    rates = "Miller"

sys_dir_unsorted = glob.glob('sys_*')
list_sys_dirs = sorted(sys_dir_unsorted, key = lambda x: int(x.split('_')[-1]))

for sys_dir in list_sys_dirs:
    simulation = CurrTempSimulation(rates=rates,analyse_Ecoul=False, source_dir=sys_dir+'/',dest_dir='analysis/'+sys_dir)

    simulation.collect_current_data()

    simulation.get_av_current()
    simulation.plot_av_current(Tlim_low=250, Tlim_high=350,Jlim_low=10000,Jlim_high=40000,errorbar=True)

    simulation.get_av_conductivity()
    simulation.plot_av_conductivity(errorbar=True)

    simulation.get_act_energy(Tlim_low=150,Tlim_high=300)

    simulation.get_conv_analysis()
    simulation.plot_conv_analysis()
