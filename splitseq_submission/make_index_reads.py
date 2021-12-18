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
    parser.add_argument('-lib_meta', dest='lib_meta',
        help='path to library metadata from google doc') # https://docs.google.com/spreadsheets/d/1r8mA9f5PMzWgSax7mSsH0r_mQGBCYnrs9SrgY6t2jKY/edit#gid=0
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
    n_ = f.count('_')
    try:
        if n_ == 4:
            tissue, age, sex, rep, read = f.split('_')
            rep = rep.split('.')[0]
        elif n_ == 3:
            tissue, age, sex, rep = f.split('_')
            rep = rep.split('.')[0]
    except: # for some reason adrenal is named weirdly
        if n_ == 4:
            tissue, age, temp, read = f.split('_')
            temp = temp.split('.')[0]
            sex = temp[0]
            rep = temp[1]
        elif n_ == 3:
            tissue, age, temp = f.split('_')
            temp = temp.split('.')[0]
            sex = temp[0]
            rep = temp[1]
    tissue = tissue.split('/')[-1]
    return tissue, age, sex, rep

def get_sample_id(f):
    s = f.split('/')[-1]
    s = s.split('.')[0]

    n_ = f.count('_')
    if n_ == 4:
        s = '_'.join(s.split('_')[:-1])
    return s

def main():
    cwd = os.path.dirname(__file__)
    args = get_args()
    d = args.dir
    long = args.long
    if args.deep:
        lib_type = 'deep'
        depth = 'deep'
    elif args.shallow:
        lib_type = 'shallow'
        depth = 'shallow'
    else:
        lib_type = None
        raise Exception('Either --deep or --shallow required')
    exclude_depth = args.exclude_depth
    if exclude_depth:
        lib_type = False
    lib_meta = args.lib_meta
    opref = args.opref

    # read metadata
    tissue_df = pd.read_csv('{}/tissue_metadata.tsv'.format(cwd), sep='\t')
    age_df = pd.read_csv('{}/age_metadata.tsv'.format(cwd), sep='\t')
    sex_df = pd.read_csv('{}/sex_metadata.tsv'.format(cwd), sep='\t')
    rep_df = pd.read_csv('{}/rep_metadata.tsv'.format(cwd), sep='\t')
    frag_df = pd.read_csv('{}/fragment_sizes.txt'.format(cwd), sep='\t',
        header=None, names=['sample', 'frag_size'])
    tissue_df.set_index('short', inplace=True)
    age_df.set_index('short', inplace=True)
    sex_df.set_index('short', inplace=True)
    rep_df.set_index('short', inplace=True)
    lib_meta = pd.read_csv(lib_meta, sep='\t')
    lib_meta['# predicted nuclei'] = lib_meta['# predicted nuclei'].astype('int')

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

        n_ = f.count('_')
        tissue, age, sex, rep = get_metadata(f)
        sample_id = get_sample_id(f)
        print(sample_id)
        # print(tissue)
        # print(age)
        # print(sex)
        # print(rep)

        # assemble metadata for this sample
        temp = tissue_df.loc[tissue]
        temp = pd.concat([temp, sex_df.loc[sex]], axis=0).to_frame()
        temp = pd.concat([temp, age_df.loc[age]], axis=0)
        temp = pd.concat([temp, rep_df.loc[int(rep)]], axis=0)
        temp = temp.transpose()

        # add biosample
        biosample = temp.copy(deep=True)
        tissue_long = temp.tissue_desc.values[0]
        if tissue != 'HC':
            desc = 'B6Cast F1 {} {} {} {}'.format(temp.age_desc.values[0], temp.model_organism_sex.values[0],
                                                  temp.tissue_desc.values[0], temp.rep_desc.values[0])
        else:
            desc = 'B6Cast F1 {} {} {} {}'.format(temp.age_desc.values[0], temp.model_organism_sex.values[0],
                                                temp.tissue_desc.values[0], temp.rep.values[0])

        biosamp_alias = 'ali-mortazavi:biosample_{}_{}_{}_{}'.format(temp.age_desc.values[0], temp.model_organism_sex.values[0],
                                                             temp.tissue_desc.values[0], temp.rep.values[0])
        biosample['date_obtained'] = lib_meta.loc[lib_meta['Sample ID'] == sample_id, 'Date shipped'].values[0]
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

        if tissue_long == 'gastrocnemius':
            tissue_long = 'gastroc'
        if lib_type:
            sample_name = tissue_long+'_'+lib_type
        else:
            sample_name = tissue_long+'_'+depth
        frag_size = frag_df.loc[frag_df['sample']==sample_name, 'frag_size'].values[0]

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
        lib['average_fragment_size'] = frag_size
        lib['construction_method'] = 'Parse Single Cell Whole Transcriptome Kit'
        lib['nucleic_acid_starting_quantity'] = lib_meta.loc[lib_meta['Sample ID'] == sample_id, '# predicted nuclei'].values[0]
        lib['nucleic_acid_starting_quantity_units'] = 'cells'

        lib['lab'] = lab
        lib['award'] = award

        cols = ['aliases', 'biosample', 'nucleic_acid_term_name',
                'documents', 'construction_method', 'nucleic_acid_starting_quantity',
                'nucleic_acid_starting_quantity_units', 'average_fragment_size',
                'strand_specificity',
                'lab', 'award']
        lib = lib[cols]

        lib_sub = pd.concat([lib_sub, lib], axis=0)

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

        cols = ['aliases', 'library', 'experiment', 'biological_replicate_number', 'technical_replicate_number']
        rep = rep[cols]

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

        file['aliases'] = file_alias
        file['dataset'] = exp_alias
        file['submitted_file_name'] = f_full_path
        file['replicate'] = rep_alias
        file['file_format'] = 'fastq'
        file['output_type'] = 'reads'
        file['lab'] = lab
        file['award'] = award

        cols = ['aliases', 'dataset', 'submitted_file_name',
                'replicate', 'file_format', 'output_type', 'platform', 'lab', 'award']
        if long == False:
            cols.append('read_length')
            cols.append('run_type')
        file = file[cols]
        file_sub = pd.concat([file_sub, file], axis=0)

        # make R2 reads
        file.rename({'aliases': 'index_of'}, axis=1, inplace=True)
        file.output_type = 'index reads'
        file.read_length = 86
        read_struct = """{"sequence_element": "UMI", "start": 1, "end": 10},"""+\
                      """{"sequence_element": "barcode", "start": 11, "end": 18},"""+\
                      """{"sequence_element": "barcode", "start": 49, "end": 56},"""+\
                      """{"sequence_element": "barcode", "start": 79, "end": 86}"""
        file['read_structure'] = read_struct
        file[['aliases_prefix', 'aliases_suffix']] = file.index_of.str.split('_', n=1, expand=True)
        file['aliases'] = file.aliases_prefix+'_r2_'+file.aliases_suffix

        if n_ == 4:
            file['file_pref'] = file.submitted_file_name.str.rsplit('.', n=2, expand=True)[0]
            file['file_pref'] = file.submitted_file_name.str.rsplit('_', n=1, expand=True)[0]
        elif n_ == 3:
            file['file_pref'] = file.submitted_file_name.str.rsplit('.', n=2, expand=True)[0]

        file.submitted_file_name = file.file_pref+'_R2.fastq.gz'

        file.drop(['aliases_prefix', 'aliases_suffix',
                   'file_pref'], axis=1, inplace=True)
        file_sub = pd.concat([file_sub, file], axis=0)

    # drop non-index reads
    file_sub = file_sub.loc[file_sub.output_type == 'index reads']

    # remove duplicates
    file_sub = file_sub.loc[~file_sub.duplicated()]

    # drop unnecessary columns
    drop = ['tissue_desc', 'age_desc', 'rep_desc', 'rep']
    biosamp_sub.drop(drop, axis=1, inplace=True)

    # for exp, also drop duplicate aliases
    exp_sub.drop_duplicates(inplace=True)

    # save each table
    if not long:
        opref = '{}_sr'.format(opref)
    else:
        opref = '{}_lr'.format(opref)

    # file
    fname = opref+'_r2_file.tsv'
    file_sub.drop('run_type', axis=1, inplace=True)
    s = file_sub.to_csv(index=False, sep='\t').replace('""', '"')
    s = s.replace('\t"', '\t')
    s = s.replace('""\t', '\t')
    s = s.replace('"\n', '\n')
    ofile = open(fname, 'w')
    ofile.write(s)
    ofile.close()

if __name__ == '__main__':
    main()
