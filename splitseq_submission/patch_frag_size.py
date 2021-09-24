import pandas as pd
import sys

# file to modify
f = sys.argv[1]

# sample name
sample = sys.argv[2]

df = pd.read_csv(f, sep='\t', usecols=[0])
df.rename({'aliases': 'record_id'}, axis=1, inplace=True)
df2 = pd.read_csv('/Users/fairliereese/Documents/programming/mortazavi_lab/bin/lab_pipelines/splitseq_submission/fragment_sizes.txt',
    sep='\t', header=None, names=['sample', 'average_fragment_size'])
n = df2.loc[df2['sample'] == sample, 'average_fragment_size'].values[0]

df['average_fragment_size'] = n
df.to_csv(f, index=False, sep='\t')
