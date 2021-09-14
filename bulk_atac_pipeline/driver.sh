#!/bin/bash
adapters=/data/users/freese/mortazavi_lab/refs/nextera_adapters/adapters.fa
bpath=/data/users/freese/mortazavi_lab/data/190626_bulk_ATAC/bin/
dpath=$1

set -e 

# input arguments
need_help=false
adapters='unassigned'
ffile='unassigned'
name='name'
genome='unassigned'
odir=$(pwd)

while getopts 'abho:f:g:i:v' flag; do
  case "${flag}" in
    f) fdir="${OPTARG}" ;;
    o) odir="${OPTARG}" ;;
    g) genome="${OPTARG}" ;;
    n) name="${OPTARG}" ;;
	i) btind=${OPTARG} ;;
    h) need_help=true ;;
	s) sample_names='unassigned' ;;
	b)
    *) error "Unexpected option ${flag}" ;;
    esac
done

#handle the help flag -h:
if [ "$need_help" = true ]
then
    echo "ATAC-seq footprint pipeline"
    echo ""
    echo "Usage: bash </path/to/atac_centipede.sh> [options]"
    echo ""
    echo "Tips:"
    echo "  Use absolute paths for all files and directories"
    echo "  Use the most recent versions of dependencies (bedtools, samtools)"
    echo ""
    echo "  Options:"
    echo "      -f <path to input fastqs>"
    echo "      -o <path to output dir>"
    echo "      -n <name> "
    echo " 		-i <path to bowtie2 index>"
    echo "      -g <reference genome in fasta format> "
    echo " 		-s <sample names as in fastq filenames, comma-separated>"
    echo "      -h <help>"
    echo ""



# # concatenate same sample/read from different lanes into the whole file
# python concatenate_lanes.py

# trim adapters
#mkdir trimmed_fastqs
#qsub trim_adapters.sh

# map to mm10
# based on doi:10.1038/nature14590 Single-cell chromatin accessibility reveals principles of regulatory variation. 
#mkdir alignments 

# sample 1
#qsub map_sample1.sh

## sample 2
#qsub map_sample2.sh

# filter out chrM reads and PCR dupes
# sample1
#bash ${bpath}filter_reads.sh /data/users/freese/mortazavi_lab/data/190626_bulk_ATAC/alignments/ATAC1_S1_unique.bam
# sample2
#bash ${bpath}filter_reads.sh /data/users/freese/mortazavi_lab/data/190626_bulk_ATAC/alignments/ATAC2_S2_unique.bam

# call peaks
mkdir peaks 
bash ${bpath}call_peaks.sh "/data/users/freese/mortazavi_lab/data/190626_bulk_ATAC/alignments/ATAC1_S1_unique_filtered.bam"

bash ${bpath}call_peaks.sh "/data/users/freese/mortazavi_lab/data/190626_bulk_ATAC/alignments/ATAC2_S2_unique_filtered.bam"
