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

# get gene names from comma separated string of genes
def get_gene_names(genes):

	if ',' in genes:
		genes = genes.split(',')
	return genes

# take the single peak index of df and split it into chrom, start, stop
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

# def coord_ind_to_coord_tuple(coord):
# 	chrom = coord.split(':')[0]
# 	start = coord.split(':')[1].split('-')[0]
# 	stop = coord.split(':')[1].split('-')[1]

# 	return (chrom, start, stop)

# def coord_tuple_to_coord_ind(chrom, start, stop):
# 	return chrom+':'+str(start)+'-'+str(stop)

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

# get the peaks that lie within the input region
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

	# append to preexisting df if not empty and we have peaks 
	# we want to save in gene_peaks
	if not peaks.empty and gene_peaks is not None:
		peaks.drop(['chrom', 'start', 'stop'], inplace=True, axis=1)
		# print('peaks')
		# print(peaks.head())
		# print(peaks.columns)
		# print(peaks.index)
		# print('gene peaks')
		# print(gene_peaks.head())
		# print(gene_peaks.columns)
		# print(gene_peaks.index)
		peaks = gene_peaks.append(peaks)
	# we have found peaks of interest but not this time
	elif peaks.empty and gene_peaks is not None: 
		peaks = gene_peaks
		print('No peaks found in input region')
	# we haven't found any peaks of interest yet
	elif peaks.empty and not gene_peaks:
		print('No peaks found in input region')
		peaks = None
	print()

	# print('Returining peaks')
	# print(peaks)
	return peaks 

# get the peak counts from counts matrix and add to our df 
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
				counts = [int(float(count)) for count in line]
				if test:
					counts = counts[:5]

				# add counts for each relevant peak
				counts_df.loc[peak_id] = [counts]

	return counts_df

# display details about reads found w/i each region
# TODO make this work with just an input region as well ie 
# remove 'gene' dependencies
def peak_deets(counts_df, any_gene=False):

	# get rid of peak indices
	num_barcodes = len(counts_df.columns)
	counts_df.reset_index(level=0, drop=True, inplace=True) 
	counts_df.reset_index(inplace=True)

	# # testing
	# print(counts_df)
	# print(counts_df.index)
	
	# group by gene name and sum to get counts ∀ cell
	counts_df = counts_df.groupby('gene').sum()

	# # testing
	# print(counts_df)
	# print(counts_df.index)

	# display details for reads in any candidate region
	if not any_gene or any_gene == 'both':
		gene_peak_deets(counts_df)

	# display details for reads in each gene's candidate region
	if any_gene or any_gene == 'both':
		print('Read details for any candidate region')
		counts_df['temp'] = 0
		counts_df = counts_df.groupby(['temp']).sum()
		counts_df.reset_index(level='temp', drop=True, inplace=True)
		disp_peak_deets(counts_df)

# details about read counts for each gene
def gene_peak_deets(df):

	# loop through each gene
	for g in df.index:
		print('Read details for %s region' % g)
		g_df = df.loc[[g]]
		disp_peak_deets(g_df)	

# details about read counts for entire list of candidate regions
def disp_peak_deets(df):
	
	num_cells = len(df.columns)

	# number of cells with reads and number of reads for each cell
	have_peak = np.count_nonzero(df.values)
	num_reads = df.values.sum()

	print('%d out of %d cells (%2f pct) have reads in region' % (have_peak, num_cells, (have_peak/num_cells)*100))
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

	# radius + whole gene body
	# chrom = 'chr'+gene['seq_region_name']
	# start = min([gene['start'], gene['end']])
	# stop = max([gene['start'], gene['end']])
	# if radius: 
	# 	temp_start = start - radius
	# 	if temp_start >= 0: start = temp_start
	# 	stop = stop + radius

	# radius + TSS
	chrom = 'chr'+gene['seq_region_name']
	start = gene['start']
	if radius:
		start = start-radius
		stop = start+radius


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

	# get the peaks 
	peaks = get_peaks(args.counts_mat)
	peaks = coord_ind_to_multiind(peaks)
	print('%d peaks' % len(peaks.index))

	# get the barcodes/cells
	barcodes = get_barcodes(args.counts_mat, args.test)
	print('%d barcodes' % len(barcodes))

	# set up some variables we'll be updating while looping
	i = 0 # index into args.genes
	gene_peaks = pd.DataFrame(columns=['peak_id', 'gene'])
	gene_peaks.set_index(['peak_id', 'gene'], inplace=True)

	# loop through each gene
	for c, sr, sp in zip(chrom, start, stop):
		print()
		if args.genes:
			print('Looking for peaks in gene: %s ' % genes[i])
		print('Region: %s %d-%d' % (c, sr, sp))

		# find peaks that lie within each region
		gene_peaks = find_peaks_in_region(c, sr, sp, peaks, gene_peaks, genes[i])
		i += 1

	# get the counts for each peak for each cell from the counts mat
	counts_df = find_peak_counts(
		args.counts_mat, gene_peaks, barcodes, args.test)

	counts_df.to_csv(make_ofile_name(args.counts_mat, 'genes'))

	# are we missing any values?
	if counts_df.isnull().sum().sum() > 0:
		print('Missing %d values ' % counts_df.isnum().sum().sum())
		print('Missing some values :/')

	# print details
	peak_deets(counts_df, 'both')
	
	# print end time
	print()
	print("--- %s seconds ---" % (time.time() - start_time))

if __name__ == '__main__':
	main()