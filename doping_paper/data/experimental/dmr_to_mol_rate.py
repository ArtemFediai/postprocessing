# convert DMR to mol_rate.
# from DMR folder csv files

import csv
import os
import pandas as pd

inp_folder = 'DMR'
out_folder = 'mol_rate'
files = os.listdir(inp_folder)

def DMR_to_mol_rate(DMR):
    return DMR / (DMR + 1)


for file_name in files:
    df = pd.read_csv(inp_folder + '/' + file_name, header=None)
    print(df.head())
    print('before:\n', df[0])
    df[0] = DMR_to_mol_rate(df[0])
    print('after:\n', df[0])

    print('WAIT')

    df.to_csv(out_folder + '/' + file_name, index=False)