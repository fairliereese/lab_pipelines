"""
Created on Sun Jul  7 07:17:14 2019

@author: alim
Adapted by Fairlie Reese
"""

import sys
import argparse
from collections import defaultdict
import numpy as np
import pandas as pd

def get_args():

	desc= 'Takes in a list of peak files to merge and outputs a merged file'
	parser = argparse.ArgumentParser(description=desc)

	parser.add_argument("-peaks", "-p", dest="pfile", 
		help = "CSV file with filepaths of peakfiles to merge")
	parser.add_argument("--o", dest="oprefix",
		help = "output merged peak file prefix, with filepath")

	args = parser.parse_args()
	return args

def get_peakfiles(pfile):

	with open(pfile, 'r') as infile:
		for line in infile:
			line = line.replace('\n', '')
			peakfiles = line.split(',')
	return peakfiles

# add peakfile to df 
def make_peakfile_df(peakfiles):

	df_list = []
	for i, pfile in enumerate(peakfiles): 
		temp = pd.read_csv(pfile, sep='\t', 
			names=['chrom', 'start', 'stop'],
			usecols=[0,1,2])
		df_list.append(temp)

	df = pd.concat(df_list, axis=0, ignore_index=True)
	
	return df

# merge peaks!
def merge_peaks(chrom_peaks):
	
	# for each chromosome
	for c, peaks in chrom_peaks.items():
		outlist = []
		# print(np.shape(peaks))
		prev_start = peaks[0][0]
		# print(prev_start)
		prev_stop = peaks[0][1]

		# loop through each peak
		for (start, stop) in peaks[1:]: 
			if prev_stop < start: 
				outlist.append((prev_start, prev_stop))
				prev_start = start
				prev_stop = stop
			elif stop <= prev_stop: 1
			else: prev_stop = stop
		# add the last peak
		outlist.append((prev_start, prev_stop))

		# replace original peaks with merged peaks in peak dict
		chrom_peaks[c] = outlist

	return chrom_peaks

# make a dictionary of peak start/end coords for each chromosome
def get_chrom_peaks(df):

	chrom_peaks = defaultdict()
	for c in df.chrom.unique():
		chrom_peaks[c] = df.loc[(df.chrom == c)][['start', 'stop']].values.tolist()
		chrom_peaks[c].sort()

	return chrom_peaks

# 
def write_peaks(merged_peaks, merge_file):

	with open(merge_file, 'w') as ofile:
		for c, peaks in merged_peaks.items():
			for p in peaks:
				ofile.write('\t'.join([str(c), str(p[0]), str(p[1])])+'\n')

#
def gen_ofile_name(oprefix):
	oprefix += '_merged.bed'
	return oprefix

def main():
	opts = get_args()
	peakfiles = get_peakfiles(opts.pfile)	

	df = make_peakfile_df(peakfiles)
	chrom_peaks = get_chrom_peaks(df)
	# print(chrom_peaks['chr1'][:10])
	chrom_peaks = merge_peaks(chrom_peaks)

	# write merged peaks to output file
	ofile = gen_ofile_name(opts.oprefix)
	write_peaks(chrom_peaks, ofile)

if __name__ == '__main__':
	main()






