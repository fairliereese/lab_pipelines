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
from collections import defaultdict
import time
import copy

def get_args():

	desc = 'Provides information about peaks in a given region from counts mat'
	parser = argparse.ArgumentParser(description=desc)

	parser.add_argument('--coords', '--c', dest='coords',
		help='coords to look for peak in format chr1:1-100')
	parser.add_argument('--genes', '--g', dest='genes',
		help='comma separated gene names to look for peaks in.')
	parser.add_argument('-matrix', '-m', dest='counts_mat',
		help='counts matrix to look in')
	parser.add_argument('--species', '--s', dest='species',
		help='species for gene name. requires --genes. default mouse.',
		default='mus musculus')
	parser.add_argument('--radius', '--r', dest='radius', type=int, 
		help='radius to look for peaks around genes. requires --genes. '+
		'default 5kb', default=5000)
	parser.add_argument("--test", dest="test", default=False, action='store_true',
		help = "only run on 5 barcodes so you can validate performance")

	args = parser.parse_args()

	# check for argument concordance
	if not args.coords and not args.genes:
		sys.exit('Missing coordinates --coords or gene name --genes')
	elif args.coords and args.genes:
		sys.exit('Please use only --coords or --genes, not both')

	# if using gene-related things, is --genes used?
	if args.species and not args.genes:
		sys.exit('--species requires --genes input')
	if args.radius and not args.genes and args.coords:
		args.radius = None
	if args.radius and not args.genes:
		sys.exit('--radius requires --genes input')

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
	# df.drop(columns='peak_ind', inplace=True)
	# df.set_index(['chrom', 'start', 'stop'], inplace=True)

	return df

def coord_ind_to_coord_tuple(coord):
	chrom = coord.split(':')[0]
	start = coord.split(':')[1].split('-')[0]
	stop = coord.split(':')[1].split('-')[1]

	return (chrom, start, stop)

def coord_tuple_to_coord_ind(chrom, start, stop):
	return chrom+':'+str(start)+'-'+str(stop)

# get the peak coordinates from the counts matrix
def get_peaks(counts_mat):

	df = pd.read_csv(counts_mat, usecols=[0])
	return df

# get the barcodes for each cell from the counts matrix
def get_barcodes(counts_mat, test):

	col_inds = pd.read_csv(counts_mat, nrows=1).columns.tolist()
	col_inds = col_inds[1:]

	# test
	if test:
		col_inds = col_inds[:5]

	return col_inds

def find_peaks_in_region(chrom, start, stop, peaks, gene_peaks, gene=None):

	# filter for relevant chromosome
	peaks = peaks.loc[peaks.chrom == chrom]

	# filter for relevant locations
	start_peaks = peaks[peaks.start.between(start, stop)]
	stop_peaks = peaks[peaks.stop.between(start, stop)]
	peaks = start_peaks.merge(stop_peaks, how='outer',
		on=['chrom', 'start', 'stop', 'peak_ind'])

	# add gene name if we're working with genes
	if gene: peaks['gene'] = gene

	# 
	# print('Peaks for gene %s ' % gene)
	# print(peaks)

	# append to preexisting df if not empty
	if not peaks.empty:
		peaks = gene_peaks.append(peaks, sort=False)
	else: 
		print('No peaks found in input region')
	print()

	return peaks 

def find_peak_counts(counts_mat, peaks, barcodes, test):

	# get peaks in list form
	peak_list = peaks['peak_ind'].values.tolist()

	# counts for each peak of interest in region
	col_ind = barcodes
	counts_df = pd.DataFrame(
		index=[peaks.peak_ind, peaks.gene], columns=col_ind)

	# read in the counts matrix one line at a time
	with open(counts_mat, 'r') as infile:
		for line in infile:
			peak_id = line.split(',')[0]

			# we found a peak we want to add
			if peak_id in peak_list:

				line = line.strip().split(',')[1:]
				# counts = [count for count in line]
				counts = [int(float(count)) for count in line]
				if test:
					counts = counts[:5]

				# add counts for each relevant peak
				# print('len counts %d '%len(counts))
				# print(type(counts))
				# print('len columns in df %d '%len(counts_df.columns))
				# print(counts_df)
				# print(peak_id)
				# print(peak_id in counts_df.index)
				# print(counts_df.index)
				# print(counts)
				counts_df.loc[peak_id] = [counts]
				# print(counts_df.head())

	return counts_df

# display details about reads found w/i each region
# TODO make this work with just an input region as well ie 
# remove 'gene' dependencies
def peak_deets(counts_df, any_gene=False):

	# mess with indices
	counts_df.reset_index(level=0, drop=True, inplace=True) 
	counts_df.reset_index(inplace=True)
	
	# group by gene name and sum to get counts ∀ cell
	counts_df = counts_df.groupby('gene').sum()

	if any_gene == 'both':

		# loop through each gene
		for g in counts_df.index:
			print('Read details for %s region' % g)
			df = counts_df.loc[g]
			print(df)
			print(type(df))
			disp_peak_deets(df)

		df = counts_df.copy(deep=True)
		df['temp'] = 0
		df = df.groupby(['temp']).sum()
		df.reset_index(level='temp', drop=True, inplace=True)
		print(type(df))
		disp_peak_deets(df)

	# count any gene into score
	elif any_gene:
		counts_df['temp'] = 0
		df = counts_df.groupby(['temp']).sum()
		df.reset_index(level='temp', drop=True, inplace=True)
		disp_peak_deets(df)

	# look at each gene individually
	elif any_gene == False:
		# loop through each gene
		for g in counts_df.index:
			print('Read details for %s region' % g)
			df = counts_df.loc[g]
			disp_peak_deets(df)



def disp_peak_deets(df):

		# number of cells
		print(type(df))
		num_barcodes = len(df.columns)

		# get number of cells that have reads and number of reads ∀ cell
		have_peak = np.count_nonzero(df.values)
		# print(have_peak)
		# have_peak = g_df.apply(lambda x: np.count_nonzero(x.values.tolist()), axis=1).values[0]
		num_reads = df.values.sum()
		# print(num_reads)
		# num_reads = g_df.apply(lambda x: sum(x.values.tolist()), axis=1).values[0]
		# print(have_peak)

		print('%d out of %d cells (%2f pct) have reads in region' % (have_peak, num_barcodes, have_peak/num_barcodes))
		print('Total of %s reads' % num_reads)
		print()

# get coordinates from input gene name 
def get_coords_from_gene_name(gene, species, radius=None):

	# species = 'mus musculus'
	try:
		gene = ensembl_rest.symbol_lookup(species=species,
		symbol=gene)
	except:
		sys.exit('Something went wrong with ENSEMBL gene name query')
	
	chrom = 'chr'+gene['seq_region_name']
	start = min([gene['start'], gene['end']])
	stop = max([gene['start'], gene['end']])

	if radius: 
		temp_start = start - radius
		if temp_start >= 0: start = temp_start
		stop = stop + radius

	return (chrom, start, stop)

def main():

	# time execution
	start_time = time.time()

	args = get_args()

	if args.genes:
		genes = get_gene_names(args.genes)
		if isinstance(genes, list):
			chroms = []
			starts = []
			stops = []
			for g in genes:
				(chrom, start, stop) = get_coords_from_gene_name(
					g, args.species, args.radius)
				chroms.append(chrom)
				starts.append(start)
				stops.append(stop)
			chrom = chroms
			start = starts
			stop = stops
		else:
			(chrom, start, stop) = get_coords_from_gene_name(
				genes, args.species, args.radius)
	else: 
		(chrom,loc) = args.coords.split(':')
		(start,stop) = loc.split('-')
		(start,stop) = (int(start),int(stop))

	# some output for the user to look at
	print()
	try:
		print('Experiment '+ args.counts_mat.split('/')[-1].split('_')[0])
	except:
		print('Experiment '+ args.counts_mat.split('_')[0])
	print('------------------------------------')


	## new way of doin it

	# get the peaks 
	peaks = get_peaks(args.counts_mat)
	peaks = coord_ind_to_multiind(peaks)
	print('%d peaks' % len(peaks.index))

	# get the barcodes
	barcodes = get_barcodes(args.counts_mat, args.test)
	print('%d barcodes' % len(barcodes))

	# loop through each gene
	i = 0
	col_inds = copy.deepcopy(barcodes)
	col_inds.append('peak_ind')
	col_inds.append('gene')
	gene_peaks = pd.DataFrame(columns=col_inds)
	for c, sr, sp in zip(chrom, start, stop):
		print()
		if args.genes:
			print('Looking for peaks in gene: %s ' % genes[i])
		print('Region: %s %d-%d' % (c, sr, sp))
		# find peaks for each region and display info
		gene_peaks = find_peaks_in_region(c, sr, sp, peaks, gene_peaks, genes[i])
		# print(gene_peaks)
		i += 1

	counts_df = find_peak_counts(
		args.counts_mat, gene_peaks, barcodes, args.test)
	counts_df.to_csv(make_ofile_name(args.counts_mat, 'genes'))
	if counts_df.isnull().sum().sum() > 0:
		print('Missing %d values ' % counts_df.isnumm().sum().sum())
		print('Missing some values :/')
	peak_deets(counts_df, 'both')
	# print end time

	print()
	print("--- %s seconds ---" % (time.time() - start_time))
	print()

	## old way of doin it
	# # read counts matrix in
	# print('Loading in counts matrix...')
	# df = pd.read_csv(args.counts_mat, sep=',')
	# df = coord_ind_to_multiind(df)

	# # some output for the user to look at
	# print()
	# try:
	# 	print('Experiment '+ args.counts_mat.split('/')[-1].split('_')[0])
	# except:
	# 	print('Experiment '+ args.counts_mat.split('_')[0])
	# print('------------------------------------')

	# # do we have more than one set of coords?
	# i = 0
	# for c, sr, sp in zip(chrom, start, stop):
	# 	print()
	# 	if args.genes:
	# 		print('Looking for peaks in gene: %s ' % genes[i])
	# 	i += 1
	# 	# find peaks and display info
	# 	peaks_df = find_peaks(c, sr, sp, df)
	# 	if peaks_df.empty: 
	# 		print('No peaks found in input region')
	# 		continue
	# 	peaks_df.to_csv(make_ofile_name(args.counts_mat, genes[i-1]))

	# 	peak_deets(peaks_df)


if __name__ == '__main__':
	main()