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

	desc = 'Takes in a bam and a barcode file and adds the barcode to the read name'
	parser = argparse.ArgumentParser(description=desc)

	parser.add_argument("-bam", "-b", dest="bamfile",
		help = "bam file from which peaks were generated")
	parser.add_argument("-barcodes", "-bc", dest="barcode_file", 
		help = "aboveKneeBarcodes file")
	parser.add_argument("-exp", "-e", dest="exp", 
		help = "Experiment ID")
	parser.add_argument("--o", dest="oprefix",
		help = "output merged peak file prefix, with filepath")

	args = parser.parse_args()
	return args

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

def main():

	# time execution
	start_time = time.time()

	args = get_args()
	exp = args.exp
	bam = args.bamfile
	barcode_file = args.barcode_file

	# get the barcodes
	barcodes = get_barcodes(barcode_file, exp)
	print(barcodes)




	# print end time
	print()
	print("--- %s seconds ---" % (time.time() - start_time))

if __name__ == '__main__':
	main()
