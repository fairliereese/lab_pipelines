import pandas as pd
import os
import sys

# usage: python get_file_ids.py <encode file submission filename>

ifile = sys.argv[1]
opref = sys.argv[2]

df = pd.read_csv(ifile, sep='\t')

df = df['aliases'].to_frame()
ofile = '{}_file_ids.tsv'.format(opref)
df.to_csv(ofile, sep='\t', index=False)
