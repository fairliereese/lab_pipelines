import argparse
import pandas as pd
import re

def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-sam', dest = 'samfile',
        help='Input sam file')
    parser.add_argument('-opref', dest='oprefix',
        help='Output prefix for files')

    args = parser.parse_args()
    return args

def main():

    # argument parsing
    args = get_args()
    samfile = args.samfile
    oprefix = args.oprefix

#     # empty DF for found mismatches
#     mm_df = pd.DataFrame(columns=['read_name', 'ref_nt', 'seq_nt', 'ind'])

    # NTs for valid mismatches
    pat = re.compile('[ACGTacgt]')

    # open output file for reading and add header
    mismatch_file = '{}_mismatches.tsv'.format(oprefix)
    outfile = open(mismatch_file, 'w')
    entry = '\t'.join(['read_name', 'ref_nt', 'seq_nt', 'ind'])+'\n'
    outfile.write(entry)

    # open and read the file
    infile = open(samfile, 'r')
    for i, line in enumerate(infile):

        if i % 400000 == 0:
            print('Processed {} reads'.format(i/4))

        # get the 4 line components of each read
        # look up the "modulo" operator if you're confused about this
        if i % 4 == 0:
            header = line.strip().split()
        elif i % 4 == 1:
            read_seq = line.strip()
        elif i % 4 == 2:
            match_seq = line.strip()
        elif i % 4 == 3:
            ref_seq = line.strip()

            # if there is a mismatch of some sort in this read and it's uniquely mapped
            if ' ' in match_seq and int(header[4]) < 256:

                # get the locations in the reference and read where the mismatches are
                inds = [m.span()[0] for m in re.finditer('[\ ]', match_seq)]

                # loop through each mismatch and get the NTs at each mismatch position
                for ind in inds:
                    ref_nt = ref_seq[ind]
                    seq_nt = read_seq[ind]

                    # if it's a true mismatch and not clipping, insertion, or deletion
                    # add it to the DF of found mismatches
                    if pat.search(ref_nt) and pat.search(seq_nt):
                        entry = '\t'.join([header[0], ref_nt, seq_nt, str(ind)])+'\n'
                        outfile.write(entry)

    infile.close()
    outfile.close()

    # count the number of mismatches per NT combination per read and
    # pivot table to format correctly
    mm_df = pd.read_csv(mismatch_file, sep='\t')
    mm_df['mismatch'] = mm_df.ref_nt+mm_df.seq_nt
    mm_df.drop(['ref_nt', 'seq_nt'], axis=1, inplace=True)
    mm_df = mm_df.groupby(by=['read_name', 'mismatch']).count().reset_index()
    mm_df.rename({'ind': 'count'}, axis=1, inplace=True)
    mm_df = mm_df.pivot(index='read_name', columns='mismatch', values='count')
    mm_df.fillna(0, inplace=True)

    ofile = '{}_mismatch_summary.tsv'.format(oprefix)
    mm_df.to_csv(ofile, sep='\t')

if __name__ == main():
    main()
