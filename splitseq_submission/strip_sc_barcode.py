#!/usr/bin/python3
from argparse import ArgumentParser
import re
import sys

def main(cmdline=None):
    parser = make_parser()
    args = parser.parse_args(cmdline)

    if args.input is None:
        instream = sys.stdin
    else:
        instream = open(args.input, 'rt')

    if args.output is None:
        outstream = sys.stdout
    else:
        outstream = open(args.output, 'wt')

    reformat_fastq(instream, outstream)

    if args.input is not None:
        instream.close()
    if args.output is not None:
        outstream.close()


def reformat_fastq(instream, outstream):
    header_re = re.compile('^@.+_VH00582')
    for line_number, line in enumerate(instream):
        if line_number % 4 == 0:
            # we are on the header line
            outstream.write(header_re.sub('@VH00582', line))
        # if one needs to alter the quality score header
        # elif (line_number + 2) % 4 == 0:
        #     outstream.write("+foo\n")
        else:
            # we are on not the sequence
            outstream.write(line)


def make_parser():
    parser = ArgumentParser()
    parser.add_argument(
        '-i',
        '--input',
        default=None,
        help="Input filename, if not provided read from standard in")
    parser.add_argument(
        "-o",
        "--output",
        default=None,
        help="Output filename if not provided write to stdout")
    return parser


if __name__ == "__main__":
    main()
