import argparse
import pandas as pd
import numpy as np

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--f", dest = "f",
        help = "read length file. one per line")
    parser.add_argument("--id", dest = "id",
        help = "id of dataset")

    args = parser.parse_args()
    return args

def compute_N50(read_lengths):
    """Given an array of read lengths, compute the N50. This metric  represents
        the read length at which half of the bases in the entire fastq are
        accounted for.
        Args:
            read_lengths(:class:`numpy.ndarray`): Array of read lengths.
        Returns:
              int: The N50 read length.
    """
    half_of_bases = read_lengths.sum()/2.0
    read_lengths.sort()
    bases_counted = 0
    for read_len in read_lengths:
        bases_counted += read_len
        if bases_counted >= half_of_bases:
            return read_len

def main():
    args = parse_args()
    df = pd.read_csv(args.f, header=None, names=['read_len'])
    read_lengths = np.array(df.read_len.tolist())
    n50 = compute_N50(read_lengths)
    print('{} N50:'.format(args.id))
    print(n50)

if __name__ == '__main__': main()
