import pandas as pd
import argparse 

def get_args():

	desc = 'Formats gene peaks csv output from check_peak.py into a '+\
		'format that upsetr can deal with'
	parser = argparse.ArgumentParser(description=desc)
	parser.add_argument('-f', '--gene_csv', dest='gene_peaks', 
		help='check_peak.py output _gene csv')
	args = parser.parse_args()
	return args

def make_ofile_name(matfile, extra):
	fname = matfile.split('.csv')[0]
	fname += '_'
	fname += extra
	fname += '.csv'
	return fname

def main():

	# get args
	args = get_args()

	# load csv
	df = pd.read_csv(args.gene_peaks)

	# groupby gene and binarize
	df = df.groupby('gene').sum()
	df[df.columns] = (df[df.columns] > 0).astype(int)

	# # remove peak_ind column
	# df.drop('peak_ind', inplace=True)

	# transpose
	df = df.transpose()

	# write to new csv
	df.to_csv(make_ofile_name(args.gene_peaks, 'upsetr'))

if __name__ == '__main__':
	main()