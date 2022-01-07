from argparse import ArgumentParser
import os
import subprocess
from collections import defaultdict
import pandas as pd
from itertools import repeat
import multiprocessing

def get_args():
    parser = ArgumentParser()
    parser.add_argument('-o', dest='opref',
        help='output directory to save split fastqs')
    parser.add_argument('-cell_meta', dest='cell_meta',
        help='path to library metadata from google doc')
    parser.add_argument('-sp_fastq', dest='sp_fastq',
        help='path to splitpipe fastq (usually named "single_cells_barcoded_head.fastq")')
    parser.add_argument('-f1', dest='fastq_1',
        help='unprocessed read 1 fastq')
    parser.add_argument('-f2', dest='fastq_2',
        help='unprocessed read 2 fastq')
    parser.add_argument('-t', dest='threads',
        help='number of threads to run on')
    parser.add_argument('--adrenal', dest='adrenal', action='store_true',
        default=False, help='whether or not to assume adrenal-style sample naming')
    args = parser.parse_args()
    return args

def create_sample_fastqs(sample, bcs, opref, fq, fq_r1, fq_r2):
    reads_oname = '{}{}_names.txt'.format(opref, sample)
    awk_str = "cat {} | grep ^@ | awk '{{if (".format(fq)
    for i, bc in enumerate(bcs):
        if i == len(bcs)-1:
            awk_str = awk_str + """substr($0,18,8)=="{}") print $0}}' | awk '{{split($0,a,"_"); print a[4]}}' | awk '{{split($0,a," "); print a[1]}}' > {}""".format(bc, reads_oname)

        else:
            awk_str += 'substr($0,18,8)=="{}"||'.format(bc)

    # run awk command and write to output file by sample
    subprocess.run(awk_str, shell=True)

    # run seqtk to find these reads
    r1_fq_oname = '{}{}.fastq'.format(opref, sample)
    r2_fq_oname = '{}{}_R2.fastq'.format(opref, sample)

    seqtk_1_str = 'seqtk subseq {} {} > {}'.format(fq_r1, reads_oname, r1_fq_oname)
    seqtk_2_str = 'seqtk subseq {} {} > {}'.format(fq_r2, reads_oname, r2_fq_oname)

    subprocess.run(seqtk_1_str, shell=True)
    subprocess.run(seqtk_2_str, shell=True)

    # delete read name files
    os.remove(reads_oname)

    # gzip fastqs
    gzip_1_str = 'gzip {}'.format(r1_fq_oname)
    gzip_2_str = 'gzip {}'.format(r2_fq_oname)

    subprocess.run(gzip_1_str, shell=True)
    subprocess.run(gzip_2_str, shell=True)

def main():
    args = get_args()
    meta_file = args.cell_meta
    opref = args.opref
    fq = args.sp_fastq
    fq_r1 = args.fastq_1
    fq_r2 = args.fastq_2
    threads = int(args.threads)
    adrenal = args.adrenal

    # get bc1 randhex / dt info
    d = os.path.dirname(__file__)
    bc_file = '{}/bc_dt_randhex.tsv'.format(d)
    bc_df = pd.read_csv(bc_file, sep='\t')

    # get sample :: well pairings
    bc_df = bc_df.melt(id_vars='well', value_vars=['bc1_dt', 'bc1_randhex'])
    bc_df.rename({'value': 'bc1'}, axis=1, inplace=True)
    bc_df.loc[bc_df.well == 0]

    # read in cell metadata
    df = pd.read_csv(meta_file)

    # limit to unique well + sample combinations
    # there should be max 48
    df = df[['sample', 'rnd1_well']].drop_duplicates()
    print('Number of unique well + sample combinations from cell_metadata: {}'.format(len(df.index)))

    # merge in sample data with bc1 data
    df = df.merge(bc_df, how='left',
                  left_on='rnd1_well',
                  right_on='well')

    # fix adrenal names
    if adrenal:
        df['sample_new'] = df['sample'].str[:-1]+'_'+df['sample'].str[-1]
        # print(df.head())
        df.drop('sample', axis=1, inplace=True)
        df.rename({'sample_new': 'sample'}, axis=1, inplace=True)
        # print(df.head())

    # create a dictionary mapping sample id :: valid bc1s
    sample_bc_dict = defaultdict(list)
    for ind, entry in df.iterrows():
        sample_bc_dict[entry['sample']].append(entry.bc1)

    # print(sample_bc_dict)
    tuples = list([(key,item) for key, item in sample_bc_dict.items()])
    samples = [i[0] for i in tuples]
    bcs = [i[1] for i in tuples]

    # # just limit to 4 samples for testing purposes
    # samples = samples[:4]
    # bcs = bcs[:4]
    # print(samples)
    # print(bcs)

    # make sure number of threads is compatible with the system
    max_cores = multiprocessing.cpu_count()
    if threads > max_cores:
    	threads = max_cores

    # write and run awk commands, then subset orig. fastqs by
    # corresponding read names
    with multiprocessing.Pool(threads) as pool:
        pool.starmap(create_sample_fastqs, zip(samples, bcs,
                     repeat(opref), repeat(fq),
                     repeat(fq_r1), repeat(fq_r2)))

if __name__ == '__main__':
    main()
