import pandas as pd
import os
import sys

# usage: python get_file_ids.py <encode file submission filename r1> <encode file submission filename r2> <opref>

ifile_r1 = sys.argv[1]
ifile_r2 = sys.argv[2]
opref = sys.argv[3]

df = pd.read_csv(ifile_r1, sep='\t')
df = df['aliases'].to_frame()

ofile = '{}_file_ids.tsv'.format(opref)
df.to_csv(ofile, sep='\t', index=False)

df = pd.read_csv(ifile_r2, sep='\t')
df = df['aliases'].to_frame()

df.to_csv(ofile, sep='\t', index=False, mode='a', header=None)


# dump path and id to new file
ofile = '{}_reupload_fastq.tsv'.format(opref)

df = pd.read_csv(ifile_r1, sep='\t')
df = df[['aliases', 'submitted_file_name']]

df.to_csv(ofile, sep='\t', index=False)

df = pd.read_csv(ifile_r2, sep='\t')
df = df[['aliases', 'submitted_file_name']]

df.to_csv(ofile, sep='\t', index=False, mode='a', header=None)
