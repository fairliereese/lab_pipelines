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
    parser.add_argument('--lib_type', dest='lib_type',
        help='additional details about exp such as "deep", "shallow", "ont_match_deep"',
        default=False)
    parser.add_argument('--deep', dest='deep', action='store_true',

    args = parser.parse_args()
    return args

def main():
    args = get_args()
    lib_meta = args.lib_meta
    folder = args.dir
    opref = args.opref
    long = args.long
    lib_type = args.lib_type

    # lib_meta = '8mo_1_shallow_metadata.tsv'
    # opref = '8mo_1_shallow'
    # long = False
    # lib_type = 'shallow'
    # folder = '/share/crsp/lab/seyedam/share/Heidi_Liz/bl6_5x_cast_ctx_hipp/fastq/shallow/'

    # get file locations
    ext = r'{}*R2.fastq.gz'.format(folder)
    locs = []
    sample_ids = []
    for f in glob.glob(ext):
        f_full_path = f
        f = os.path.basename(f)
        locs.append(f_full_path)
        sample_ids.append(f.split('_R2.fastq.gz')[0])

    file_loc = pd.DataFrame(data=sample_ids,
                             columns=['sample_id'])
    file_loc['r2_loc'] = locs
    file_loc['r1_loc'] = file_loc['r2_loc'].str.split('_R2.fastq.gz', expand=True)[0]+'.fastq.gz'

    # metadata
    meta_dir = os.path.dirname(__file__)
    tissue = pd.read_csv('{}/tissue_metadata.tsv'.format(meta_dir), sep='\t')
    tissue.drop('donor', axis=1, inplace=True)
    geno = pd.read_csv('{}/genotype_metadata.tsv'.format(meta_dir), sep='\t')
    sex = pd.read_csv('{}/sex_metadata.tsv'.format(meta_dir), sep='\t')
    age = pd.read_csv('{}/age_metadata.tsv'.format(meta_dir), sep='\t')
    rep = pd.read_csv('{}/rep_metadata.tsv'.format(meta_dir), sep='\t')

    # biosample
    df = pd.read_csv(lib_meta, sep='\t')
    df['tissue'] = df['Sample ID'].str.split('_', expand=True)[1]
    df['genotype'] = df['Sample ID'].str.split('_', expand=True)[0]
    df['replicate'] = df['Sample ID'].str.split('_', expand=True)[2].astype(int)
    df = df.merge(tissue, how='left', left_on='tissue',
        right_on='short')
    df = df.merge(geno, how='left', left_on='genotype',
        right_on='short')
    df = df.merge(sex, how='left', left_on='Sex',
        right_on='short')
    df = df.merge(age, how='left', left_on='Timepoint',
        right_on='short')
    df = df.merge(rep, how='left', left_on='replicate',
        right_on='short')
    df['date_obtained'] = df['Date shipped']
    df['description'] = df.desc_genotype+' '+\
                        df.age_desc+' '+\
                        df.model_organism_sex+' '+\
                        df.tissue_desc+' '+\
                        df.rep.astype(str)
    df['aliases'] = 'ali-mortazavi:biosample_'+\
                        df.alias_genotype+'_'+\
                        df.age_desc+'_'+\
                        df.model_organism_sex+'_'+\
                        df.tissue_desc+'_'+\
                        df.rep.astype(str)
    df['lab'] = 'ali-mortazavi'
    df['award'] = 'UM1HG009443'

    cols = ['aliases', 'alias_genotype', 'age_desc',
            'model_organism_sex', 'tissue_desc', 'rep']
    biosamp_merge = df[cols].copy(deep=True)
    biosamp_merge.rename({'aliases': 'biosample'}, axis=1, inplace=True)

    cols = ['biosample_ontology', 'organism', 'subcellular_fraction_term_name',
            'donor', 'source', 'model_organism_sex', 'model_organism_age',
            'model_organism_age_units', 'date_obtained', 'description',
            'aliases', 'lab', 'award']
    biosamp = df[cols]
    fname = '{}_biosample.tsv'.format(opref)
    biosamp.to_csv(fname, sep='\t', index=False)

    # experiment
    df['assay_term_name'] = 'single-cell RNA sequencing assay'
    if not long:
        df['description'] = 'Short-read Split-seq '+\
                            df.desc_genotype+' '+\
                            df.age_desc+' '+\
                            df.model_organism_sex+' '+\
                            df.tissue_desc
        df['aliases'] = 'ali-mortazavi:exp_sr_'+\
                            df.alias_genotype+'_'+\
                            df.age_desc+'_'+\
                            df.model_organism_sex+'_'+\
                            df.tissue_desc
    else:
        df['description'] = 'Long-read Split-seq '+\
                            df.desc_genotype+' '+\
                            df.age_desc+' '+\
                            df.model_organism_sex+' '+\
                            df.tissue_desc
        df['aliases'] = 'ali-mortazavi:exp_lr_'+\
                            df.alias_genotype+'_'+\
                            df.age_desc+'_'+\
                            df.model_organism_sex+'_'+\
                            df.tissue_desc

    cols = ['aliases', 'alias_genotype', 'age_desc',
            'model_organism_sex', 'tissue_desc']
    exp_merge = df[cols].copy(deep=True)
    exp_merge.rename({'aliases': 'experiment'}, axis=1, inplace=True)

    cols = ['aliases', 'biosample_ontology', 'description',\
            'assay_term_name', 'lab', 'award']
    exp = df[cols].copy(deep=True)

    if lib_type:
        exp['description'] = exp.description+'_'+lib_type
        exp['aliases'] = exp.aliases+'_'+lib_type
    exp.drop_duplicates(inplace=True)

    fname = '{}_experiment.tsv'.format(opref)
    exp.to_csv(fname, sep='\t', index=False)


    # library
    if not long:
        df['aliases'] = 'ali-mortazavi:lib_sr_'+\
                            df.alias_genotype+'_'+\
                            df.age_desc+'_'+\
                            df.model_organism_sex+'_'+\
                            df.tissue_desc+'_'+\
                            df.rep.astype(str)
    else:
        df['aliases'] = 'ali-mortazavi:lib_lr_'+\
                            df.alias_genotype+'_'+\
                            df.age_desc+'_'+\
                            df.model_organism_sex+'_'+\
                            df.tissue_desc+'_'+\
                            df.rep.astype(str)

    df = df.merge(biosamp_merge, how='left', on=['alias_genotype', 'age_desc',
            'model_organism_sex', 'tissue_desc', 'rep'])

    if lib_type:
        df['description'] = df.description+'_'+lib_type
        df['aliases'] = df.aliases+'_'+lib_type
        if 'shallow' in lib_type:
            df['nucleic_acid_starting_quantity'] = df['# predicted nuclei for shallow datasets']
        elif 'deep' in lib_type:
            df['nucleic_acid_starting_quantity'] = df['# predicted nuclei for deep datasets']

    df['nucleic_acid_term_name'] = 'RNA'
    df['documents'] = 'ali-mortazavi:Split-seq_computational_protocol_v1.0,ali-mortazavi:split-seq-v1'
    df['construction_method'] = 'Parse Single Cell Whole Transcriptome Kit'
    df['nucleic_acid_starting_quantity_units'] = 'cells'
    df['strand_specificity'] = 'unstranded'
    df['average_fragment_size'] = df['Fragment size (bp)']

    cols = ['aliases', 'alias_genotype', 'age_desc',
            'model_organism_sex', 'tissue_desc', 'rep']
    lib_merge = df[cols].copy(deep=True)
    lib_merge.rename({'aliases': 'library'}, axis=1, inplace=True)

    cols = ['aliases', 'biosample', 'nucleic_acid_term_name', \
            'documents', 'construction_method', 'nucleic_acid_starting_quantity', \
            'nucleic_acid_starting_quantity_units', 'strand_specificity', \
            'average_fragment_size', 'lab', 'award']
    lib = df[cols].copy(deep=True)
    fname = '{}_library.tsv'.format(opref)
    lib.to_csv(fname, sep='\t', index=False)

    # replicate
    if long == False:
        df['aliases'] = 'ali-mortazavi:rep_sr_'+\
                      df.alias_genotype+'_'+\
                      df.age_desc+'_'+\
                      df.model_organism_sex+'_'+\
                      df.tissue_desc+'_'+\
                      df.rep.astype(str)
    else:
        df['aliases'] = 'ali-mortazavi:rep_lr_'+\
                      df.alias_genotype+'_'+\
                      df.age_desc+'_'+\
                      df.model_organism_sex+'_'+\
                      df.tissue_desc+'_'+\
                      df.rep.astype(str)

    if lib_type:
        df['aliases'] + df.aliases+'_{}'.format(lib_type)

    df = df.merge(exp_merge, how='left', on=['alias_genotype', 'age_desc',
            'model_organism_sex', 'tissue_desc'])
    df = df.merge(lib_merge, how='left', on=['alias_genotype', 'age_desc',
            'model_organism_sex', 'tissue_desc', 'rep'])
    df['biological_replicate_number'] = df['rep']
    df['technical_replicate_number'] = 1

    cols = ['aliases', 'alias_genotype', 'age_desc',
            'model_organism_sex', 'tissue_desc', 'rep']
    rep_merge = df[cols].copy(deep=True)
    rep_merge.drop_duplicates(inplace=True)
    rep_merge.rename({'aliases': 'replicate'}, axis=1, inplace=True)

    cols = ['aliases', 'library', 'experiment',
            'biological_replicate_number', 'technical_replicate_number']
    rep = df[cols]
    rep.drop_duplicates(inplace=True)
    fname = '{}_rep.tsv'.format(opref)
    rep.to_csv(fname, sep='\t', index=False)

    # r1 fastqs
    if long == False:
        df['aliases'] = 'ali-mortazavi:fastq_r1_sr_'+\
                      df.alias_genotype+'_'+\
                      df.age_desc+'_'+\
                      df.model_organism_sex+'_'+\
                      df.tissue_desc+'_'+\
                      df.rep.astype(str)

        df['read_length'] = 115
        df['run_type'] = 'single-ended'
        df['platform'] = 'encode:NextSeq2000'
    else:
        df['aliases'] = 'ali-mortazavi:fastq_lr_'+\
                      df.alias_genotype+'_'+\
                      df.age_desc+'_'+\
                      df.model_organism_sex+'_'+\
                      df.tissue_desc+'_'+\
                      df.rep.astype(str)
        df['platform'] = 'encode:PacBio_sequel_II'

    if lib_type:
        df['aliases'] + df.aliases+'_{}'.format(lib_type)

    df.drop('replicate', axis=1, inplace=True)
    df = df.merge(rep_merge, how='left', on=['alias_genotype', 'age_desc',
            'model_organism_sex', 'tissue_desc', 'rep'])
    df['dataset'] = df['experiment']
    df['file_format'] = 'fastq'
    df['output_type'] = 'reads'

    cols = ['aliases', 'alias_genotype', 'age_desc',
            'model_organism_sex', 'tissue_desc', 'rep']
    r1_merge = df[cols].copy(deep=True)
    r1_merge.rename({'aliases': 'index_of'}, axis=1, inplace=True)

    df = df.merge(file_loc, how='left', left_on='Sample ID',
            right_on='sample_id')
    df['submitted_file_name'] = df['r1_loc']

    cols = ['aliases', 'dataset', 'submitted_file_name',
            'replicate', 'file_format', 'output_type', 'platform', 'lab', 'award']
    if long == False:
        cols.append('read_length')
        cols.append('run_type')
    file = df[cols].copy(deep=True)
    file.drop_duplicates(inplace=True)
    fname = '{}_r1_fastq.tsv'.format(opref)
    file.to_csv(fname, sep='\t', index=False)

    # r2 fastqs
    if long == False:
        df['aliases'] = 'ali-mortazavi:fastq_r2_sr_'+\
                      df.alias_genotype+'_'+\
                      df.age_desc+'_'+\
                      df.model_organism_sex+'_'+\
                      df.tissue_desc+'_'+\
                      df.rep.astype(str)

        df['read_length'] = 86
        df['run_type'] = 'single-ended'
        df['platform'] = 'encode:NextSeq2000'

        if lib_type:
            df['aliases'] + df.aliases+'_{}'.format(lib_type)

        df['file_format'] = 'fastq'
        df['output_type'] = 'index reads'
        df['read_structure'] = '{"sequence_element": "UMI", "start": 1, "end": 10},{"sequence_element": "cell barcode", "start": 11, "end": 18},{"sequence_element": "cell barcode", "start": 49, "end": 56},{"sequence_element": "cell barcode", "start": 79, "end": 86}'
        cols =['alias_genotype', 'age_desc',
        'model_organism_sex', 'tissue_desc', 'rep']
        df = df.merge(r1_merge, how='left', on=['alias_genotype', 'age_desc',
                'model_organism_sex', 'tissue_desc', 'rep'])
        df['submitted_file_name'] = df['r2_loc']

        cols = ['aliases', 'dataset', 'submitted_file_name',
                'replicate', 'file_format', 'output_type',
                'platform', 'lab', 'award', 'read_length',
                'run_type', 'index_of', 'read_structure']
        file = df[cols].copy(deep=True)
        file.drop_duplicates(inplace=True)
        fname = '{}_r2_fastq.tsv'.format(opref)
        file.to_csv(fname, sep='\t', index=False)









# if __name__ == '__main__':
#     main()
