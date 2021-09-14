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
		help = "only run on 200 peaks so you can validate performance")
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

# convert bam file to sam file so it can be parsed
def bam_to_sam(bamfile, oprefix):
	samfile = oprefix+'.sam'
	pysam.view('-h', bamfile, '-o', samfile, catch_stdout=False)
	return samfile

# read peakfile into dataframe
def read_peakfile(pfile, test):
	df = pd.read_csv(pfile, sep='\t', 
		names=['chrom', 'start', 'stop', 'x', 'y'])

	df.drop(['x', 'y'], axis=1, inplace=True)
	# df['peak_ind'] = df.apply(lambda x: x.chrom+':'+str(x.start)+'-'+str(x.stop), axis=1)
	# df.set_index('peak_ind', inplace=True)
	# print(df.head())

	# testing
	if test:
		df = df.head(200)

	print('Using %d peaks' % len(df.index))

	return df

# read barcodes into dataframe
def get_barcodes(barcode_file, exp):
	barcodes = []

	with open(barcode_file, 'r') as infile:
		for i,line in enumerate(infile):
			if i != 0:
				bc = line.strip().split(',')[0]
				bc = bc.split('tagged_')[1]
				# df['_'.join([exp, bc])] = []
				barcodes.append(bc)

	exps = [exp for i in range(len(barcodes))]
	ind = ['_'.join([exp, b]) for b in barcodes]
	df = pd.DataFrame(index=ind)
	df['barcode'] = pd.Series(barcodes).values
	df['exp'] = pd.Series(exps).values

	return df

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
	df = pd.DataFrame(data, index=row_inds, columns=col_inds)

	# # add combined coordinates field to use as index after populating
	# peak_ind = [chroms[i]+':'+str(starts[i])+'-'+str(stops[i])
	# 	for i in range(len(chroms))]
	# df['peak_ind'] = peak_ind

	return df

# count the number of reads in each peak
def get_peak_counts(peak_df, bfile, count_df, exp, args):

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

	# drop chrom, start, stop cols to save space, reindex with peak_ind
	count_df.reset_index(inplace=True)
	count_df.rename({'level_0': 'chrom', 'level_1': 'start', 'level_2': 'stop'},
		inplace=True, axis=1)
	# add peak index
	count_df['peak_ind'] = count_df.apply(
		lambda x: x.chrom+':'+str(x.start)+'-'+str(x.stop), axis=1)
	count_df.drop(['chrom', 'start', 'stop'], inplace=True, axis=1)
	count_df.set_index('peak_ind', inplace=True)

	# make this a sparse df
	if args.sparse_csv or args.sparse_anndata:
		print('Converting to sparse matrix...')
		count_df = count_df.astype(pd.SparseDtype(float), copy=False)	


	print(count_df.head())
	print(count_df.index)
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
	barcodes = get_barcodes(args.barcode_file, exp)

	# initialize anndata matrix of 0s 
	count_df = init_count_df(barcodes.barcode.values.tolist(), peak_df, exp)
	# print(count_df.head())
	count_df = get_peak_counts(peak_df, bfile, count_df, exp, args)

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



	# print(count_df.head())


	# count_df.to_csv(exp+'_counts.csv')




if __name__ == '__main__':
	main()

