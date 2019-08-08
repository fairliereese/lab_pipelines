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

	# parser.add_argument('--coords', '--c', dest='coords',
	# 	help='coords to look for peak in format chr1:1-100')
	# parser.add_argument('--genes', '--g', dest='genes',
	# 	help='comma separated gene names to look for peaks in.')
	parser.add_argument('-gc', '--genes_coords', dest='genes',
		help='multiline file of comma separated list of genes and/or '+
		'coordinates to look for. optional <label=cell_type> can be '+
		'be put at the end of line. end format will look like '+
		'gene1,gene2,label=myoblast<newline>gene3,gene4,label=myotube')
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
	parser.add_argument("--oprefix", '--o', dest='oprefix', type=str,
		help='prefix to add to output file. do not put path in this',
		default=None)

	args = parser.parse_args()

	return args

def make_ofile_name(matfile, gene, prefix=None):
	fname = matfile.split('.csv')[0]
	if prefix:
		fname += '_'
		fname += prefix
	fname += '_'
	fname += gene
	fname += '.csv'
	return fname

# get gene names from comma separated string of genes
def get_gene_names(genes):

	if ',' in genes:
		genes = genes.split(',')

	if isinstance(genes, list):
		return genes
	else:
		return [genes]

# given an input file of format
# gene1,gene2,gene3,label=MT
# or
# chr1:123-456,chr2:789-101
# parse and extract genes to list and labelled dictionary
def get_gene_names_file(fname): 

	genes_df = pd.DataFrame(columns=['label', 'gene'])
	with open(fname, 'r') as infile:
		for i, line in enumerate(infile):

			# get the group label for each input group of genes
			# and the genes
			try:
				genes = line.strip().split(',')
				group_label = line.strip().split('label=')[-1]
				genes = genes[:-1]
			except:
				group_label = i

			# add this line's genes to df 
			g_df = pd.DataFrame(columns=['label', 'gene'])
			g_df['gene'] = genes
			g_df['label'] = group_label
			genes_df = genes_df.append(g_df)


			# # add to dictionary
			# labelled_genes[group_label]	= genes	

	# return labelled_genes, genes
	return genes_df

# take the single peak index of df and split it into chrom, start, stop
def coord_ind_to_multiind(df):

	df['chrom'] = df.apply(
		lambda x: x.peak_ind.split(':')[0], axis=1)
	df['start'] = df.apply(
		lambda x: int(x.peak_ind.split(':')[1].split('-')[0]), axis=1)
	df['stop'] = df.apply(
		lambda x: int(x.peak_ind.split(':')[1].split('-')[1]), axis=1)

	return df

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

# display details about reads found w/i each region
# TODO make this work with just an input region as well ie 
# remove 'gene' dependencies
def peak_deets(counts_df):


	# testing
	# print(counts_df)
	# print(counts_df.index)

	# get rid of peak indices
	num_barcodes = len(counts_df.columns)
	counts_df.reset_index(level='peak_ind', drop=True, inplace=True) 
	counts_df.reset_index(inplace=True)

	# # testing
	# print(counts_df)
	# print(counts_df.index)
	
	# group by gene name and sum to get counts ∀ cell
	counts_df = counts_df.groupby(['gene', 'label']).sum()
	# counts_df.reset_index(level='label')

	# # testing
	# print(counts_df)
	# print(counts_df.index)

	# deets_df
	deets_df = pd.DataFrame(columns=['label', 'num_reads', 'percent_positive'])

	# display details for reads in each gene's candidate region
	temp = counts_df.reset_index(level='label', drop=True)
	deets_df = gene_peak_deets(temp, deets_df)


	# print(counts_df)

	# display details for reads in each labelled set of genes
	temp = counts_df.groupby(level='label').sum()
	temp.reset_index(level='label', inplace=True)
	temp.sort_values('label', inplace=True)
	# print(temp)
	# print(temp)
	for label in temp.label.unique():
		print('Read details for label {0} '.format(label))
		label_df = temp.loc[temp.label == label]
		label_df.set_index('label', inplace=True)
		# print(label_df)
		deets_df = disp_peak_deets(label_df, label, deets_df)

	return deets_df

# details about read counts for each gene
def gene_peak_deets(df, deets_df):

	# loop through each gene
	for g in df.index:
		print('Read details for %s region' % g)
		g_df = df.loc[[g]]
		deets_df = disp_peak_deets(g_df, g, deets_df)	

	return deets_df

# details about read counts for entire list of candidate regions
def disp_peak_deets(df, label, deets_df):

	# print(df)
	# print(df.values)
	
	num_cells = len(df.columns)

	# number of cells with reads and number of reads for each cell
	have_peak = np.count_nonzero(df.values)
	num_reads = df.values.sum()

	print('%d out of %d cells (%2f pct) have reads in region' % (have_peak, num_cells, (have_peak/num_cells)*100))
	print('Total of %s reads' % num_reads)
	print()

	temp = pd.DataFrame(data=[[label, num_reads, (have_peak/num_cells)*100]],
		columns=['label', 'num_reads', 'percent_positive'])

	deets_df = deets_df.append(temp)
	return deets_df




# # get coordinates from input gene name 
# def get_coords_from_gene_name(gene, species, radius=None):

# 	# species = 'mus musculus'
# 	try:
# 		gene = ensembl_rest.symbol_lookup(species=species,
# 		symbol=gene)
# 	except:
# 		sys.exit('Something went wrong with ENSEMBL gene name query')

# 	# radius + whole gene body
# 	chrom = 'chr'+gene['seq_region_name']
# 	start = min([gene['start'], gene['end']])
# 	stop = max([gene['start'], gene['end']])
# 	if radius: 
# 		temp_start = start - radius
# 		if temp_start >= 0: start = temp_start
# 		stop = stop + radius

# 	# # radius + TSS
# 	# chrom = 'chr'+gene['seq_region_name']
# 	# start = gene['start']
# 	# if radius:
# 	# 	stop = start+radiuss
# 	# 	start = start-radius


	# return (chrom, start, stop)

def get_coords(x, args):

	g = x.gene 
	if '-' in g and ':' in g: 
			(chrom,loc) = g.split(':')
			(start,stop) = loc.split('-')
			(start,stop) = (int(start),int(stop))
	# gene name
	else:
		try:
			gene = ensembl_rest.symbol_lookup(species=args.species,
				symbol=g)
		except:
			sys.exit('Something went wrong with ENSEMBL gene name query')

		# radius + TSS
		chrom = 'chr'+gene['seq_region_name']
		temp1 = gene['start'] - 5000 
		temp2 = gene['start'] + 5000
		start = min(temp1, temp2)
		stop = max(temp1, temp2)

		# # testing
		# print(g)
		# print('start: {}'.format(gene['start']))
		# print('region start: {}'.format(start))
		# print('region end: {}'.format(stop))

		# # radius + whole gene body
		# chrom = 'chr'+gene['seq_region_name']
		# start = min([gene['start'], gene['end']])
		# stop = max([gene['start'], gene['end']])
		# if args.radius: 
		# 	temp_start = start - args.radius
		# 	if temp_start >= 0: start = temp_start
		# 	stop = stop + args.radius

		# # testing
		# print('new start: {}'.format(start))
		# print('new stop: {}'.format(stop))


	return pd.Series([chrom, start, stop])
 
def get_coords_df(genes_df, args):

	new_cols = genes_df.apply(get_coords, args=(args,), axis=1)
	new_cols.columns = ['chrom', 'start', 'stop']

	# make a new dataframe 
	genes_df['chrom'] = new_cols['chrom']
	genes_df['start'] = new_cols['start']
	genes_df['stop'] = new_cols['stop']

	return genes_df

def find_relevant_peaks_df(x, peaks):

	peak_ind = x.chrom+':'+str(x.start)+'-'+str(x.stop)

	print('Finding peaks for {0} in region {1}'.format(x.gene, peak_ind))

	region_peaks = pd.DataFrame(columns=['label', 'gene', 'peak_ind'])

	# filter for relevant chromosome
	peaks = peaks.loc[peaks.chrom == x.chrom]

	# filter for locations
	start_peaks = peaks[peaks.start.between(x.start, x.stop)]
	stop_peaks = peaks[peaks.stop.between(x.start, x.stop)]
	peaks = start_peaks.merge(stop_peaks)

	# add gene and label
	peaks['gene'] = x.gene
	peaks['label'] = x.label

	# remove chrom, start, stop
	peaks.drop(['chrom', 'start', 'stop'], inplace=True, axis=1)

	# are there even any peaks?
	if peaks.empty:
		print('    No peaks found in this input region')
		return peaks

	# print(peaks)

	return peaks

def find_relevant_peaks(genes_df, peaks):

	# init df
	region_peaks = pd.DataFrame(columns=['label', 'gene', 'peak_ind'])
	region_peaks.set_index(['label', 'gene', 'peak_ind'], inplace=True)

	# get peak ids of relevant peaks
	temp = genes_df.apply(find_relevant_peaks_df, args=(peaks,), axis=1)

	# concatenate all the series returned from apply
	for i, s in temp.iteritems():
		region_peaks = region_peaks.append(s)

	return region_peaks

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
def get_peak_counts(counts_mat, region_peaks, barcodes, test):

	# get peaks in list form
	peak_list = region_peaks['peak_ind'].values.tolist()

	# counts for each peak of interest in region
	col_ind = barcodes
	counts_df = pd.DataFrame(
		index=[region_peaks.peak_ind, region_peaks.gene, region_peaks.label], columns=col_ind)

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

def main():

	# time execution
	start_time = time.time()

	# input args
	args = get_args()
	genes = args.genes
	counts_mat = args.counts_mat
	oprefix = args.oprefix

	# some output for the user to look at
	print()
	try:
		print('Experiment '+ counts_mat.split('/')[-1].split('_')[0])
	except:
		print('Experiment '+ counts_mat.split('_')[0])
	print('------------------------------------')

	# get genes/coords
	genes_df = get_gene_names_file(genes)
	genes_df = get_coords_df(genes_df, args)
	# print(genes_df)

	# get the peaks 
	peaks = get_peaks(counts_mat)
	# print(peaks)
	peaks = coord_ind_to_multiind(peaks)
	print('%d peaks' % len(peaks.index))

	# get the barcodes/cells
	barcodes = get_barcodes(counts_mat, args.test)
	print('%d barcodes' % len(barcodes))

	# get the peaks associated with each target region
	region_peaks = find_relevant_peaks(genes_df, peaks)

	# get the counts for each relevant region
	counts_df = get_peak_counts(
		counts_mat, region_peaks, barcodes, args.test)
	# print(counts_df)

	# save counts info 
	counts_df.to_csv(
		make_ofile_name(counts_mat, 'region_peaks', oprefix))

	# display details
	print()
	summary_df = peak_deets(counts_df)
	summary_df.to_csv(
		make_ofile_name(counts_mat, 'region_summary', oprefix), index=False)

	# print end time
	print()
	print("--- %s seconds ---" % (time.time() - start_time))

if __name__ == '__main__':
	main()