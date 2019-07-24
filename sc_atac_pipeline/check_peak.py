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
	parser.add_argument('--genes', '--g', dest='genes',
		help='comma separated gene names to look for peaks in')
	parser.add_argument('-matrix', '-m', dest='counts_mat',
		help='counts matrix to look in')

	args = parser.parse_args()
	return args

def make_ofile_name(matfile, gene):
	fname = matfile.split('.csv')[0]
	fname += '_'
	fname += gene
	fname += '.csv'
	return fname

def get_gene_names(genes):

	if ',' in genes:
		genes = genes.split(',')
	return genes

def coord_ind_to_multiind(df):

	df['chrom'] = df.apply(
		lambda x: x.peak_ind.split(':')[0], axis=1)
	df['start'] = df.apply(
		lambda x: int(x.peak_ind.split(':')[1].split('-')[0]), axis=1)
	df['stop'] = df.apply(
		lambda x: int(x.peak_ind.split(':')[1].split('-')[1]), axis=1)
	df.drop(columns='peak_ind', inplace=True)
	# df.set_index(['chrom', 'start', 'stop'], inplace=True)

	return df

# finds peaks within a region defined by the user
def find_peaks(chrom, start, stop, df):

	print('Searching for peaks in region: %s %d-%d' % (chrom, start, stop))
	coord_range = [i for i in range(start, stop+1)]

	peaks = df.loc[df.chrom == chrom]
	print('After filtering on chromosome number...')
	print(peaks.chrom.unique())
	peaks = peaks[peaks.apply(
		lambda x: x.start in coord_range or x.stop in coord_range, axis=1)]
	print('After filtering on peak location')
	print(peaks.chrom.unique())

	peaks = peaks.copy()
	
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
	if not args.coords and not args.genes:
		sys.exit('Missing coordinates --coords or gene name --genes')
	elif args.coords and args.genes:
		sys.exit('Please use only --coords or --genes, not both')
	elif args.genes:
		genes = get_gene_names(args.genes)
		if isinstance(genes, list):
			chroms = []
			starts = []
			stops = []
			for g in genes:
				(chrom, start, stop) = get_coords_from_gene_name(g)
				chroms.append(chrom)
				starts.append(start)
				stops.append(stop)
			chrom = chroms
			start = starts
			stop = stops
		else:
			(chrom, start, stop) = get_coords_from_gene_name(genes)
	else: 
		(chrom,loc) = args.coords.split(':')
		(start,stop) = loc.split('-')
		(start,stop) = (int(start),int(stop))

	# read counts matrix in
	print('Loading in counts matrix...')
	df = pd.read_csv(args.counts_mat, sep=',')
	df = coord_ind_to_multiind(df)

	# some output for the user to look at
	print()
	try:
		print('Experiment '+ args.counts_mat.split('/')[-1].split('_')[0])
	except:
		print('Experiment '+ args.counts_mat.split('_')[0])
	print('------------------------------------')

	# do we have more than one set of coords?
	i = 0
	for c, sr, sp in zip(chrom, start, stop):
		print()
		if args.genes:
			print('Looking for peaks in gene: %s ' % genes[i])
		i += 1
		# find peaks and display info
		peaks_df = find_peaks(c, sr, sp, df)
		if peaks_df.empty: 
			print('No peaks found in input region')
			continue
		peaks_df.to_csv(make_ofile_name(args.counts_mat, genes[i-1]))

		peak_deets(peaks_df)


if __name__ == '__main__':
	main()