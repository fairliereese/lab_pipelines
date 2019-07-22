import sys
import glob
import os

# make combined fastas folder
if not os.path.isdir('combined_fastqs'):
	os.mkdir('combined_fastqs')

outputs = [['ATAC1', 'S1', 'R1'], \
		   ['ATAC1', 'S1', 'R2'], \
		   ['ATAC2', 'S2', 'R1'], \
		   ['ATAC2', 'S2', 'R2']]

# make an output fasta file for each different sample combo
out_fastas = []
for o in outputs:
	out_fastas.append(open('combined_fastqs/'+'_'.join(o)+'.fastq', 'w'))

# open up each fasta and dump the contents to the big one
s = glob.glob('*fastq')
#s_split = [i.split('_')[2] for i in s]
#print(s_split)
#print(s)
#print(type(s))
s = sorted(s, key=lambda x: x.split('_')[2])
#print(s)
for file in s:
	print(file)

	# which sample file does this need to go to 
	for i, o in enumerate(outputs):
		if all(substr in file for substr in o):
			out_fasta = out_fastas[i]

	print(out_fasta.name)
	print()

	# actually open subject fasta
	with open(file, 'r') as f: 
		seqs = f.read()
		out_fasta.write(seqs)

