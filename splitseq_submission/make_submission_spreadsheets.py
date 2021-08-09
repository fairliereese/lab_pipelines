import pandas as pd
from argparse import ArgumentParser
import os
import glob

def get_args():
    parser = ArgumentParser()
    parser.add_argument('-d', dest='dir',
        help='directory where split fastqs are')
    args = parser.parse_args()

def get_metadata(f):
    tissue, age, sex, rep = f.split('_')
    rep = rep.split('.')[0]
    return tissue, age, sex, rep

def add_biosamp(tissue, age, sex, rep):

def __main__():
    args = get_args()

    # iterate through the input files
    ext = r'{}*.fastq.gz'.format(args.dir)
    for f in glob.glob(ext):

        tissue, age, sex, rep = get_metadata(f)

        biosamp = add_biosamp(tissue, age, sex, rep)

if __name__ == '__main__':
    main()
