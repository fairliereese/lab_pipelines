# adapted from laksshman (lsundaram on github)
import h5py
import pandas as pd 
import argparse


def get_args():
	desc = 'Takes a .snap file and converts it to a .csv file'
	parser = argparse.ArgumentParser(description=desc)

	parser.add_argument('-s', '--snap', dest='snapfile',
		help='Snapfile to convert from .snap to .csv')
	parser.add_argument('-bs', '--binsize', dest='bin',
		help='Binsize in base pairs to bin your data on')

	args = parser.parse_args()
	return args

def bin_snapfile(f, binsize):

	dfa = pd.DataFrame()
	try:
		dfa['binchr'] = f['AM'][str(binsize)]['binChrom'][0:]
		dfa['binstart'] = f['AM'][str(binsize)]['binStart'][0:]

		dfb = pd.DataFrame()
		dfb['idx'] = f['AM'][str(binsize)]['idx']
		dfb['idy'] = f['AM'][str(binsize)]['idy']
		dfb['count'] = f['AM'][str(binsize)]['count']

		return dfa, dfb

	except:
		print('Bins of this size do not exist yet for this snapfile')
		exit()

def get_snapfile(snapfile):
	f = h5py.File(snapfile, 'r')
	return f

def main():
	args = get_args()
	snapfile = args.snapfile
	binsize = args.bin

	f = get_snapfile()

	dfa, dfb = bin_snapfile(f, binsize)

if __name__ == '__main__': main()
