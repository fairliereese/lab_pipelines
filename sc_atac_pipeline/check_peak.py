#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jul  7 19:34:04 2019

@author: alim
Adapted by Fairlie Reese
"""
import argparse
import pandas as pd
import numpy as np
import ensembl_rest
import sys

def get_args():

	desc = 'Provides information about peaks in a given region from counts mat'
	parser = argparse.ArgumentParser(description=desc)

	parser.add_argument('--coords', '--c', dest='coords',
		help='coords to look for peak in format chr1:1-100')
	parser.add_argument('--gene_name', '--g', dest='gene_name',
		help='gene name to look for peaks in')
	parser.add_argument('-matrix', '-m', dest='counts_mat',
		help='counts matrix to look in')

	args = parser.parse_args()
	return args

# finds peaks within a region defined by the user
def find_peaks(chrom, start, stop, df):

	print('Searching for peaks in region: %s %d-%d' % (chrom, start, stop))
	coord_range = [i for i in range(start, stop+1)]

	peaks = df.loc[df.chrom == chrom]
	peaks = df[df.apply(
		lambda x: x.start in coord_range or x.stop in coord_range, axis=1)]
	
	return peaks

# display details of peaks found in region to user
def peak_deets(df):

	df.set_index(['chrom', 'start', 'stop'], inplace=True)

	df['thisisdumb'] = 0
	num_barcodes = len(df.columns)

	df = df.groupby(['thisisdumb']).sum()
	df.reset_index(level='thisisdumb', drop=True, inplace=True)

	have_peak = df.apply(lambda x: np.count_nonzero(x.values.tolist()), axis=1).values[0]
	num_reads = df.apply(lambda x: sum(x.values.tolist()), axis=1).values[0]
	# print(have_peak)

	print('%d out of %d cells (%.2f pct) have reads in region' % (have_peak, num_barcodes, have_peak/num_barcodes))
	print('Total of %d reads' % num_reads)

# get coordinates from input gene name 
def get_coords_from_gene_name(gene):

	species = 'mus musculus'
	try:
		gene = ensembl_rest.symbol_lookup(species=species,
		symbol=gene)
	except:
		sys.exit('Something went wrong with ENSEMBL gene name query')
	
	chrom = 'chr'+gene['seq_region_name']
	start = min([gene['start'], gene['end']])
	stop = max([gene['start'], gene['end']])

	return (chrom, start, stop)

def main():
	args = get_args()

	# did the user give a gene name or coords?
	if not args.coords and not args.gene_name:
		sys.exit('Missing coordinates --coords or gene name --gene_name')
	elif args.coords and args.gene_name:
		sys.exit('Please use only --coords or --gene_name, not both')
	elif args.gene_name:
		(chrom, start, stop) = get_coords_from_gene_name(args.gene_name)
	else: 
		(chrom,loc) = args.coords.split(':')
		(start,stop) = loc.split('-')
		(start,stop) = (int(start),int(stop))

	# read counts matrix in
	df = pd.read_csv(args.counts_mat, sep='\t')

	# some output for the user to look at
	print()
	try:
		print('Experiment '+ args.counts_mat.split('/')[-1].split('_')[0])
	except:
		print('Experiment '+ args.counts_mat.split('_')[0])
	print('------------------------------------')
	print()

	# find peaks and display info
	df = find_peaks(chrom, start, stop, df)
	if df.empty: 
		print('No peaks found in input region')
		return

	peak_deets(df)

if __name__ == '__main__':
	main()