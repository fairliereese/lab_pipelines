#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jul  6 21:44:29 2019

@author: alim
Adapted by Fairlie Reese
"""

import sys
import argparse
from collections import defaultdict
import numpy as np
import pandas as pd
import pysam
from multiprocessing import Process

def get_args():

	desc = 'Takes in a list of peak files to merge and outputs a merged file'
	parser = argparse.ArgumentParser(description=desc)

	parser.add_argument("--peaks", "-p", dest="peakfile", 
		help = "peak file in bed format")
	parser.add_argument("--bam", "-b", dest="bamfile",
		help = "bam file from which peaks were generated")
	parser.add_argument("--barcodes", "-bc", dest="barcode_file", 
		help = "aboveKneeBarcodes file")
	parser.add_argument("--exp", "-e", dest="exp", 
		help = "Experiment ID")
	parser.add_argument("-o", dest="oprefix",
		help = "output merged peak file prefix, with filepath")

	args = parser.parse_args()
	return args

# convert bam file to sam file so it can be parsed
def bam_to_sam(bamfile, oprefix):
	samfile = oprefix+'.sam'
	pysam.view('-h', bamfile, '-o', samfile, catch_stdout=False)
	return samfile

# read peakfile into dataframe
def read_peakfile(pfile, exp):
	df = pd.read_csv(pfile, sep='\t', 
		names=['chrom', 'start', 'stop', 'x', 'y'])

	df.drop(['x', 'y'], axis=1, inplace=True)
	df['exp'] = exp

	print('Using %d peaks' % len(df.index))

	return df

# get a list of each barcode from the barcode file
def get_barcodes(barcode_file):
	barcodes = []

	with open(barcode_file, 'r') as infile:
		for i,line in enumerate(infile):
			if i != 0:
				bc = line.strip().split(',')[0]
				bc = bc.split('tagged_')[1]
				barcodes.append(bc)

	# print(len(barcodes))
	# barcodes = list(set(barcodes))
	# print(len(barcodes))
	return barcodes

# initialize a df of 0s with a spot for each combination of peaks/barcodes/experiment
def init_count_df(barcodes, peak_df, exp):
	
	# get lists of values to index df on 
	row_inds = [peak_df.chrom.values.tolist(), 
				peak_df.start.values.tolist(),
				peak_df.stop.values.tolist()]
	col_inds = [exp+'_'+b for b in barcodes]
	# col_inds.append(barcodes)

	# make matrix of zeros 
	data = np.zeros(shape=(len(row_inds[0]), len(col_inds)))

	# init zeros df
	df = pd.DataFrame(data, index=row_inds, columns=col_inds)
	return df

# count the number of reads in each peak
def get_peak_counts(peak_df, bfile, count_df, exp):

	# some things to keep track of during iteration
	with_barcode = 0
	without_barcode = 0 
	read_counter = 0

	samfile = pysam.AlignmentFile(bfile, "rb")
	for i, peak in peak_df.iterrows():

		start = int(peak.start)
		stop = int(peak.stop)
		chrom = peak.chrom

		for read in samfile.fetch(chrom, start-1, stop):

			read_counter += 1
			if read_counter % 1000000 == 0:
				print('Processed %d reads' % read_counter)

			# get the barcode from the read
			try:
				b = read.get_tag('RG').split('_')
				b = '_'.join([b[1],b[2]])
				with_barcode += 1
				# print(b)
			except: 
				without_barcode += 1
				continue

		# update count in df
		col_ind = '_'.join([exp, b])
		count_df.loc[(chrom, start, stop), col_ind] += 1

	print('%d reads have have barcode' % with_barcode)
	print('%d reads missing barcode' % without_barcode)

	return count_df

def main():
	args = get_args()
	pfile = args.peakfile
	bfile = args.bamfile
	exp = args.exp

	# get peaks from bed file
	peak_df = read_peakfile(pfile, exp)

	# break up into each of the 23 chroms

	barcodes = get_barcodes(args.barcode_file)

	count_df = init_count_df(barcodes, peak_df, exp)
	print(count_df.head())
	count_df = get_peak_counts(peak_df, bfile, count_df, exp)
	print(count_df.head())


	count_df.to_csv(exp+'_counts.csv')




if __name__ == '__main__':
	main()

