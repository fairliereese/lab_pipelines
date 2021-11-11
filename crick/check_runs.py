import pandas as pd

df1 = pd.read_csv('210816_transfer_runs.tsv', sep='\t',
    header=None)
df2 = pd.read_csv('210816_transfer_runs_1.tsv', sep='\t',
    header=None)
df3 = pd.read_csv('211104_transfer_runs.tsv', sep='\t',
    header=None)

df1.columns = df2.columns = df3.columns = ['run_id']
df1['df1'] = True
df2['df2'] = True
df3['df3'] = True

# df1['run_id'] = df1['run_id'].astype(object)
# df2['run_id'] = df2['run_id'].astype(object)
# df3['run_id'] = df3['run_id'].astype(object)

df1 = df1.merge(df2, how='outer', on='run_id')
df1 = df1.merge(df3, how='outer', on='run_id')
