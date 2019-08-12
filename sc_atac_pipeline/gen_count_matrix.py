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
import anndata as ad
import scipy 
import time

def get_args():

	desc = 'Takes in a list of peak files to merge and outputs a merged file'
	parser = argparse.ArgumentParser(description=desc)

	parser.add_argument("-peaks", "-p", dest="peakfile", 
		help = "peak file in bed format")
	parser.add_argument("-bam", "-b", dest="bamfile",
		help = "bam file from which peaks were generated")
	parser.add_argument("-barcodes", "-bc", dest="barcode_file", 
		help = "aboveKneeBarcodes file")
	parser.add_argument("-exp", "-e", dest="exp", 
		help = "Experiment ID")
	parser.add_argument("--o", dest="oprefix",
		help = "output merged peak file prefix, with filepath")
	parser.add_argument("--test", dest="test", default=False, action='store_true',
		help = "only run on 200 peaks/5 barcodes so you can validate performance")
	parser.add_argument("--full_csv", dest="full_csv", default=False, action='store_true',
		help = "output to full csv")
	parser.add_argument("--sparse_csv", dest="sparse_csv", default=False, action='store_true',
		help = "output to sparse csv")
	parser.add_argument("--sparse_anndata", dest="sparse_anndata", default=False, action='store_true',
		help = "output to sparse anndata")
	parser.add_argument("--hd5", dest="hd5", default=False, action='store_true',
		help = "output to hd5 file")


	args = parser.parse_args()
	return args

def gen_fname(oprefix, exp, test):
	if oprefix[-1] == '/': 
		fname = oprefix+exp
	# elif '/' in oprefix: 
	# 	fname = oprefix+'/'+exp
	else: 
		fname = oprefix+exp

	if test:
		fname += '_test'
	return fname

# read peakfile into dataframe
def read_peakfile(pfile, test):
	df = pd.read_csv(pfile, sep='\t', 
		names=['chrom', 'start', 'stop', 'x', 'y'])

	df.drop(['x', 'y'], axis=1, inplace=True)

	# testing
	if test:
		df = df.head(200)

	print('Using %d peaks' % len(df.index))

	return df

# read barcodes into dataframe
def get_barcodes(barcode_file, exp, test):
	barcodes = []

	with open(barcode_file, 'r') as infile:
		for i,line in enumerate(infile):
			if i != 0:
				bc = line.strip().split(',')[0]
				bc = bc.split('tagged_')[1]
				# df['_'.join([exp, bc])] = []
				barcodes.append(bc)

	if test:
		exps = [exp for i in range(len(barcodes[:5]))]
		ind = ['_'.join([exp, b]) for b in barcodes[:5]]
		barcodes = barcodes[:5]
	else: 
		exps = [exp for i in range(len(barcodes))]
		ind = ['_'.join([exp,b]) for b in barcodes]

	df = pd.DataFrame(index=ind)
	df['barcode'] = pd.Series(barcodes).values
	df['exp'] = pd.Series(exps).values

	print('Using %d barcodes' % len(barcodes))

	return df, ind

# initialize a df of 0s with a spot for each combination of peaks/barcodes/experiment
def init_count_df(barcodes, peak_df, exp):
	
	# get lists of values to index df on 
	chroms = peak_df.chrom.values.tolist()
	starts = peak_df.start.values.tolist()
	stops = peak_df.stop.values.tolist()

	row_inds = [chroms, starts, stops]
	col_inds = [exp+'_'+b for b in barcodes]

	# make matrix of zeros 
	data = np.zeros(shape=(len(row_inds[0]), len(col_inds)))

	# init zeros df
	# df = pd.DataFrame(data, index=row_inds, columns=col_inds)
	# row-wise addition
	df = pd.DataFrame(columns=col_inds)


	return df

# count the number of reads in each peak
def get_peak_counts(peak_df, bfile, count_df, exp, args, cols):

	# some things to keep track of during iteration
	with_barcode = 0
	without_barcode = 0 
	read_counter = 0

	samfile = pysam.AlignmentFile(bfile, "rb")
	for i, peak in peak_df.iterrows():

		start = int(peak.start)
		stop = int(peak.stop)
		chrom = peak.chrom

		# set up dataframe to hold each peak's counts
		cols = np.asarray(cols)
		data = np.zeros(shape=(1,len(cols)))
		peak_row = pd.DataFrame(data, columns=cols)
		peak_ind = chrom+':'+str(start)+'-'+str(stop)
		peak_row['peak_ind'] = peak_ind
		peak_row.set_index('peak_ind', inplace=True)

		for read in samfile.fetch(chrom, start-1, stop):

			read_counter += 1
			if read_counter % 1000000 == 0:
				print('Processed %d reads' % read_counter)

			# get the barcode from the read
			try:
				b = read.get_tag('RG').split('_')
				b = '_'.join([b[1],b[2]])
				with_barcode += 1
			except: 
				without_barcode += 1
				continue

			# add updated row to df
			col_ind = '_'.join([exp, b])
			if args.test: 
				if col_ind in cols:
					peak_row.loc[:, col_ind] += 1
			else: peak_row.loc[:, col_ind] += 1
		count_df = count_df.append(peak_row)

	print('%d reads have have barcode' % with_barcode)
	print('%d reads missing barcode' % without_barcode)

	# make this a sparse df
	if args.sparse_csv or args.sparse_anndata:
		print('Converting to sparse matrix...')
		count_df = count_df.astype(pd.SparseDtype(float), copy=False)	


	print(count_df.head())
	print()

	return count_df

def make_anndata(count_df, peaks, cells):

	peaks['peak_ind'] = peaks.apply(lambda x: x.chrom+':'+str(x.start)+'-'+str(x.stop), axis=1)
	peaks.set_index('peak_ind', inplace=True)

	print(count_df.index)
	print(peaks.index)
	print(cells.index)
	df = ad.AnnData(X=count_df, obs=peaks, var=cells)
	return df

def main():

	# time execution
	start_time = time.time()

	args = get_args()
	pfile = args.peakfile
	bfile = args.bamfile
	exp = args.exp
	test = args.test

	if args.oprefix:
		fname = gen_fname(args.oprefix, exp, test)

	# get peaks from bed file
	peak_df = read_peakfile(pfile, test)

	# get barcodes from barcode file
	barcodes, col_ind = get_barcodes(args.barcode_file, exp, test)

	# initialize anndata matrix of 0s 
	count_df = init_count_df(barcodes.barcode.values.tolist(), peak_df, exp)
	count_df = get_peak_counts(peak_df, bfile, count_df, exp, args, col_ind)

	# write to csv full matrix
	if args.full_csv:
		print()
		print('Pandas full version')
		print('Writing to output file...')
		count_df.to_csv(fname+'_full_counts.csv')

	# write to csv sparse matrix
	if args.sparse_csv:
		print()
		print('Pandas sparse version')
		print('Writing to output file...')
		count_df.to_csv(fname+'_sparse_counts.csv')

	# write pandas hd5
	if args.hd5:
		print()
		print('Pandas hd5 version')
		print('Writing to output file...')
		count_df.to_hdf(fname+'_counts.hdf5', key='count_df')

	# write anndata hd5
	if args.sparse_anndata:
	# get anndata version
		print('Anndata version')
		count_df = make_anndata(count_df, peak_df, barcodes)
		print()
		print('Writing to output file...')
		count_df.write_h5ad(fname+'_sparse_counts.h5ad')

	# print end time
	print()
	print("--- %s seconds ---" % (time.time() - start_time))

if __name__ == '__main__':
	main()

