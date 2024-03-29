import pandas as pd
from argparse import ArgumentParser
import os
import glob
import sys

def get_args():
    parser = ArgumentParser()
    parser.add_argument('-d', dest='dir',
        help='directory where split fastqs are')
    parser.add_argument('-o', dest='opref',
        help='output directory to save tsvs')
    parser.add_argument('--long', dest='long', action='store_true',
        help='long read experiment')
    parser.add_argument('--deep', dest='deep', action='store_true',
        help='deep sequencing experiment')
    parser.add_argument('--shallow', dest='shallow', action='store_true',
        help='shallow sequencing experiment')
    parser.add_argument('--exclude_depth', dest='exclude_depth', action='store_true')

    args = parser.parse_args()
    return args

def get_metadata(f):
    try:
        tissue, age, sex, rep = f.split('_')
        rep = rep.split('.')[0]
    except: # for some reason adrenal is named weirdly
        tissue, age, temp = f.split('_')
        temp = temp.split('.')[0]
        sex = temp[0]
        rep = temp[1]
    tissue = tissue.split('/')[-1]
    return tissue, age, sex, rep

def get_sample_id(f):
    s = f.split('/')[-1]
    s = s.split('.')[0]
    return s

def main():
    cwd = os.path.dirname(__file__)
    args = get_args()
    d = args.dir
    long = args.long
    if args.deep:
        lib_type = 'deep'
    elif args.shallow:
        lib_type = 'shallow'
    else:
        lib_type = None
        raise Exception('Either --deep or --shallow required')
    exclude_depth = args.exclude_depth
    if exclude_depth:
        lib_type = False
    opref = args.opref

    # read metadata
    tissue_df = pd.read_csv('{}/tissue_metadata.tsv'.format(cwd), sep='\t')
    age_df = pd.read_csv('{}/age_metadata.tsv'.format(cwd), sep='\t')
    sex_df = pd.read_csv('{}/sex_metadata.tsv'.format(cwd), sep='\t')
    rep_df = pd.read_csv('{}/rep_metadata.tsv'.format(cwd), sep='\t')
    swap_df = pd.read_csv('{}/sample_swaps.tsv'.format(cwd), sep='\t')
    tissue_df.set_index('short', inplace=True)
    age_df.set_index('short', inplace=True)
    sex_df.set_index('short', inplace=True)
    rep_df.set_index('short', inplace=True)

    # globals
    lab = 'ali-mortazavi'
    award = 'UM1HG009443'

    ext = r'{}*.fastq.gz'.format(d)

    # create submission tables
    biosamp_sub = pd.DataFrame()
    exp_sub = pd.DataFrame() # only one for each pair of reps!
    lib_sub = pd.DataFrame()
    rep_sub = pd.DataFrame()
    file_sub = pd.DataFrame()

    for f in glob.glob(ext):
        f_full_path = f
        f = os.path.basename(f)

        tissue, age, sex, rep = get_metadata(f)
        sample_id = get_sample_id(f)
        print(sample_id)

        # assemble metadata for this sample
        temp = tissue_df.loc[tissue]
        temp = pd.concat([temp, sex_df.loc[sex]], axis=0).to_frame()
        temp = pd.concat([temp, age_df.loc[age]], axis=0)
        temp = pd.concat([temp, rep_df.loc[int(rep)]], axis=0)
        temp = temp.transpose()

        # add biosample
        biosample = temp.copy(deep=True)
        desc = 'B6Cast F1 {} {} {} {}'.format(temp.age_desc.values[0], temp.model_organism_sex.values[0],
                                              temp.tissue_desc.values[0], temp.rep_desc.values[0])
        biosamp_alias = 'ali-mortazavi:biosample_{}_{}_{}_{}'.format(temp.age_desc.values[0], temp.model_organism_sex.values[0],
                                                             temp.tissue_desc.values[0], temp.rep.values[0])
        biosample['description'] = desc
        biosample['aliases'] = biosamp_alias
        biosample['lab'] = lab
        biosample['award'] = award

        biosamp_sub = pd.concat([biosamp_sub, biosample], axis=0)

        # add experiment
        exp = temp.copy(deep=True)
        if long == False:
            desc = 'Short-read Split-seq B6Cast F1 {} {} {}'.format(temp.age_desc.values[0], temp.model_organism_sex.values[0],
                                                  temp.tissue_desc.values[0])
            exp_alias = 'ali-mortazavi:exp_sr_{}_{}_{}'.format(temp.age_desc.values[0], temp.model_organism_sex.values[0],
                                                         temp.tissue_desc.values[0])
        else:
            desc = 'Long-read Split-seq B6Cast F1 {} {} {}'.format(temp.age_desc.values[0], temp.model_organism_sex.values[0],
                                                  temp.tissue_desc.values[0])
            exp_alias = 'ali-mortazavi:exp_lr_{}_{}_{}'.format(temp.age_desc.values[0], temp.model_organism_sex.values[0],
                                                 temp.tissue_desc.values[0])

        if lib_type:
            desc += ' {}'.format(lib_type)
            exp_alias += '_{}'.format(lib_type)

        exp['assay_term_name'] = 'single-cell RNA sequencing assay'
        exp['description'] = desc
        exp['aliases'] = exp_alias
        exp['lab'] = lab
        exp['award'] = award

        cols = ['aliases', 'biosample_ontology', 'description',
                'assay_term_name', 'lab', 'award']
        exp = exp[cols]

        exp_sub = pd.concat([exp_sub, exp], axis=0)

        # add library
        lib = temp.copy(deep=True)

        if long == False:
            lib_alias = 'ali-mortazavi:library_sr_{}_{}_{}_{}'.format(temp.age_desc.values[0], temp.model_organism_sex.values[0],
                                                                 temp.tissue_desc.values[0], temp.rep.values[0])
            lib['documents'] = 'ali-mortazavi:Split-seq_computational_protocol_v1.0,ali-mortazavi:split-seq-v1'
            lib['strand_specificity'] = 'unstranded'
        else:
            lib_alias = 'ali-mortazavi:library_lr_{}_{}_{}_{}'.format(temp.age_desc.values[0], temp.model_organism_sex.values[0],
                                                         temp.tissue_desc.values[0], temp.rep.values[0])
            lib['documents'] = 'ali-mortazavi:LR-Split-seq_computational_protocol_v1.0,ali-mortazavi:split-seq-v1,ali-mortazavi:pacbio-split-seq-v1'
            lib['strand_specificity'] = 'forward'

        if lib_type:
            lib_alias += '_{}'.format(lib_type)

        lib['aliases'] = lib_alias
        lib['biosample'] = biosamp_alias
        lib['nucleic_acid_term_name'] = 'RNA'
        lib['construction_method'] = 'Parse Single Cell Whole Transcriptome Kit'
        lib['nucleic_acid_starting_quantity_units'] = 'cells'

        lib['lab'] = lab
        lib['award'] = award

        # add replicate
        rep = temp.copy(deep=True)

        if long == False:
            rep_alias = 'ali-mortazavi:rep_sr_{}_{}_{}_{}'.format(temp.age_desc.values[0], temp.model_organism_sex.values[0],
                                                             temp.tissue_desc.values[0], temp.rep.values[0])
        else:
            rep_alias = 'ali-mortazavi:rep_lr_{}_{}_{}_{}'.format(temp.age_desc.values[0], temp.model_organism_sex.values[0],
                                                             temp.tissue_desc.values[0], temp.rep.values[0])

        if lib_type:
            rep_alias += '_{}'.format(lib_type)

        rep['aliases'] = rep_alias
        rep['library'] = lib_alias
        rep['experiment'] = exp_alias
        rep['biological_replicate_number'] = temp.rep.values[0]
        rep['technical_replicate_number'] = 1

        rep_sub = pd.concat([rep_sub, rep], axis=0)

        # add file
        file = temp.copy(deep=True)

        if long == False:
            file_alias = 'ali-mortazavi:fastq_sr_{}_{}_{}_{}'.format(temp.age_desc.values[0], temp.model_organism_sex.values[0],
                                                             temp.tissue_desc.values[0], temp.rep.values[0])
            file['read_length'] = 115
            file['run_type'] = 'single-ended'
            file['platform'] = 'encode:NextSeq2000'

        else:
            file_alias = 'ali-mortazavi:fastq_lr_{}_{}_{}_{}'.format(temp.age_desc.values[0], temp.model_organism_sex.values[0],
                                                             temp.tissue_desc.values[0], temp.rep.values[0])
            file['platform'] = 'encode:PacBio_sequel_II'

        if lib_type:
            file_alias += '_{}'.format(lib_type)

        file['record_id'] = file_alias
        file['dataset'] = exp_alias
        file['replicate'] = rep_alias
        file['sample_id'] = sample_id

        cols = ['record_id', 'dataset', 'replicate', 'sample_id']
        file = file[cols]
        file_sub = pd.concat([file_sub, file], axis=0)

    # save each table
    if not long:
        opref = '{}_sr'.format(opref)
    else:
        opref = '{}_lr'.format(opref)

    # fix sample swaps
    swaps = pd.read_csv('{}/sample_swaps.tsv'.format(cwd), sep='\t')
    # old_records = df.merge(swaps['sample_id_old'], how='inner', left_on='sample_id', right_on='sample_id_old')
    print(swaps.columns)
    print(file_sub.columns)
    swaps = file_sub.merge(swaps[['sample_id_old', 'sample_id_new']], how='inner', left_on='sample_id', right_on='sample_id_old')
    swaps.drop('sample_id_old', axis=1, inplace=True)

    # new_records = df.merge(swaps['sample_id_new'], how='inner', left_on='sample_id_old', right_on='sample_id_new')
    # new_records.rename({'record_id_old': 'record_id_new',
    #                     'dataset'})
    swaps = swaps.merge(file_sub, how='inner', left_on='sample_id_new', right_on='sample_id', suffixes=('_new', '_old'))
    swaps.record_id_new = swaps.record_id_new+'_temp'
    swaps.rename({'record_id_new': 'aliases',
                  'dataset_new': 'dataset',
                  'replicate_new': 'replicate',
                  'record_id_old': 'record_id'}, axis=1, inplace=True)
    swaps = swaps[['record_id', 'aliases', 'dataset', 'replicate']]
    swaps.head()

    fname = '{}_file_patch_1.tsv'.format(opref)
    swaps.to_csv(fname, index=False, sep='\t')

    swaps = swaps[['aliases']]
    swaps.rename({'aliases': 'record_id'}, axis=1, inplace=True)
    swaps.to_csv('test1.csv')
    print(swaps.head())
    swaps['aliases'] = swaps.record_id.str.rsplit(pat='_', n=1, expand=True)[0]
    print(swaps.head())
    swaps.to_csv('test2.csv')
    fname = '{}_file_patch_2.tsv'.format(opref)
    swaps.to_csv(fname, index=False, sep='\t')

if __name__ == '__main__':
    main()
