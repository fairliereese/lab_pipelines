import pandas as pd
import glob
import os
import numpy as np

def get_metadata(f,
                 batch='normal'):
    """
    Return tissue, age, sex, replicate, and genotype metadata
    from a filename.

    Parameters:
        f (str): Base filename
        eight_mo (str): Whether this is 8mo or normal exp.
    """
    if batch != '8mo':
        # R2s have one extra set of characters
        if 'R2' in f:
            tissue, age, sex, rep, _ = f.split('_')
            rep = rep.split('.')[0]
        else:
            tissue, age, sex, rep = f.split('_')
            rep = rep.split('.')[0]
        geno = 'B'
    else:
        age = '8mo'
        # R2s have one extra set of characters
        if 'R2' in f:
            geno, tissue, rep, _ = f.split('_')
            rep = rep.split('.')[0]
        else:
            geno, tissue, rep = f.split('_')
            rep = rep.split('.')[0]
        sex = np.nan
    return tissue, age, sex, rep, geno

def load_ref_data():
    refs = dict()

    refs['age'] = pd.read_csv('refs/age_metadata.tsv', sep='\t')

    refs['sex'] = pd.read_csv('refs/sex_metadata.tsv', sep='\t')

    refs['tissue'] = pd.read_csv('refs/tissue_metadata.tsv', sep='\t')
    refs['tissue'].drop('donor', axis=1, inplace=True)

    refs['rep'] = pd.read_csv('refs/rep_metadata.tsv', sep='\t')
    refs['rep'].short = refs['rep'].short.astype(str)

    refs['geno'] = pd.read_csv('refs/genotype_metadata.tsv', sep='\t')

    refs['assay'] = pd.read_csv('refs/assay_metadata.tsv', sep='\t')

    # limit library metadata based on what we're looking for?
    refs['lib'] = pd.read_csv('refs/lib_metadata.tsv', sep='\t')
    cols = ['Sample ID', 'Tissue', 'depth2', 'batch', 'run_number', 'Date shipped',
            'Fragment size (bp)', '# of input cells',
            'Sex', 'Path to R1 fastq', 'assay']
    refs['lib'] = refs['lib'][cols]
    refs['lib']['Path to R1 fastq'] = refs['lib']['Path to R1 fastq'].astype(str)
    refs['lib']['dir'] = refs['lib'].apply(lambda x: os.path.dirname(x['Path to R1 fastq'])+'/', axis=1)
    refs['lib'].rename({'depth2':'depth',
                        'Date shipped': 'date_obtained',
                        'Fragment size (bp)': 'average_fragment_size',
                        '# of input cells': 'nucleic_acid_starting_quantity'}, axis=1, inplace=True)
#     refs['lib'].date_obtained = refs['lib'].date_obtained.dt.strftime('%Y-%m-%d')
    refs['lib'].date_obtained = refs['lib'].date_obtained.astype('datetime64')

    return refs

refs = load_ref_data()

temp = refs['lib'][['assay', 'batch', 'depth', 'dir', 'Sample ID']].groupby(['assay', 'batch', 'depth', 'dir']).count().reset_index()

first_lib = True
for ind, entry in temp.iterrows():
    if entry.dir == '/':
        continue

    assay = entry.assay
    batch = entry.batch
    depth = entry.depth
    d = entry.dir

    print(assay)
    print(batch)
    print(depth)
    print(d)
    print()

    # find all the fastqs
    files = []
    ext = r'{}*.fastq.gz'.format(d)
    for f in glob.glob(ext):
        files.append(f)
    files = pd.DataFrame(data=files, columns=['fname'])

    df = files.copy(deep=True)

    # add library information
    df['depth'] = depth
    df['batch'] = batch
    df['assay'] = assay

    # add consistent info
    df['lab'] = 'ali-mortazavi'
    df['award'] = 'UM1HG009443'
    df['nucleic_acid_term_name'] = 'RNA'
    df['construction_method'] = 'Parse Single Cell Whole Transcriptome Kit'
    df['nucleic_acid_starting_quantity_units'] = 'cells'
    df['technical_replicate_number'] = 1

    df['short_fname'] = df.apply(lambda x: os.path.basename(x.fname), axis=1)
    df[['tissue', 'age', 'sex', 'rep', 'geno']] = df.apply(lambda x: get_metadata(x.short_fname, batch=batch), axis=1, result_type='expand')

    if batch == 'normal' or batch == 'ont_match':
        df['Sample ID'] = df.tissue+'_'+df.age+'_'+df.sex+'_'+df.rep
        df['biological_replicate_number'] = df.rep
    elif batch == '8mo':
        df['Sample ID'] = df.geno+'_'+df.tissue+'_'+df.rep

    # tissue-related metadata
    df = df.merge(refs['tissue'], how='left', left_on='tissue', right_on='short')

    # library-related metadata
    df['tissue_merge'] = df.tissue_desc.str.capitalize()
    df = df.merge(refs['lib'], how='left', on=['Sample ID', 'depth', 'batch', 'assay'])

    # 8mo data needs to pull sex from the spreadsheet
    if batch == '8mo':
        df['sex'] = df['Sex']

    # age-related metadata
    df = df.merge(refs['age'], how='left', left_on='age', right_on='short')

    # sex-related metadata
    df = df.merge(refs['sex'], how='left', left_on='sex', right_on='short')

    # replicate-related metadata
    df = df.merge(refs['rep'], how='left', left_on='rep', right_on='short')

    # 8mo data needs to pull rep from metadata
    if batch == '8mo':
        df['biological_replicate_number'] = df.replicate
        df['rep'] = df['replicate'].astype('str')

    # genotype-related metadata
    df = df.merge(refs['geno'], how='left', left_on='geno', right_on='short')

    # assay-related metadata
    print(df.columns)
    df = df.merge(refs['assay'], how='left', left_on='assay', right_on='short')

    # biosample
    b = df.copy(deep=True)

    # create the alias
    b['aliases'] = 'ali-mortazavi:biosamp_'+b.alias_genotype+ \
                   '_'+b.tissue_desc+ \
                   '_'+b.age_desc+ \
                   '_'+b.model_organism_sex+ \
                   '_'+b.rep

    # create the description
    b[['desc_genotype', 'tissue_desc', 'age_desc', 'model_organism_sex', 'rep_desc']]
    b['description'] = b.desc_genotype+ \
                       ' '+b.tissue_desc+ \
                       ' '+b.age_desc+ \
                       ' '+b.model_organism_sex+ \
                       ' '+b.rep

    cols = ['biosample_ontology', 'organism', 'subcellular_fraction_term_name',
            'donor', 'source', 'model_organism_sex', 'model_organism_age',
            'model_organism_age_units', 'date_obtained', 'description',
            'aliases', 'lab', 'award', 'Sample ID']
    b = b[cols]
    b.drop_duplicates(inplace=True)

    # merge in alias with orig. df to have access for libraries
    b_alias = b[['aliases', 'Sample ID']]
    b_alias.rename({'aliases': 'biosample'}, axis=1, inplace=True)
    df = df.merge(b_alias, how='left', on='Sample ID')
    b.drop('Sample ID', axis=1, inplace=True)
    b.drop_duplicates(inplace=True)

    fname = 'biosample.tsv'
    if first_lib:
        b.to_csv(fname, sep='\t', index=False)
    else:
        b.to_csv(fname, sep='\t', index=False, header=None, mode='a')


    # experiment
    e = df.copy(deep=True)

    # create the alias
    e['aliases'] = 'ali-mortazavi:experiment_'+e.alias_assay+ \
                   '_'+e.alias_genotype+ \
                   '_'+e.tissue_desc+ \
                   '_'+e.age_desc+ \
                   '_'+e.model_organism_sex+ \
                   '_'+e.batch+ \
                   '_'+e.depth

    # create the description
    e['description'] = e.desc_assay+ \
                       ' '+e.desc_genotype+ \
                       ' '+e.tissue_desc+ \
                       ' '+e.age_desc+ \
                       ' '+e.model_organism_sex+ \
                       ' '+e.batch+ \
                       ' '+e.depth

    cols = ['aliases', 'biosample_ontology', 'description', \
            'assay_term_name', 'lab', 'award', 'Sample ID']
    e = e[cols]
    e.drop_duplicates(inplace=True)

    # merge in alias with orig. df to have access for libraries
    e_alias = e[['aliases', 'Sample ID']]
    e_alias.rename({'aliases': 'experiment'}, axis=1, inplace=True)
    df = df.merge(e_alias, how='left', on='Sample ID')
    e.drop('Sample ID', axis=1, inplace=True)
    e.drop_duplicates(inplace=True)

    fname = 'experiment.tsv'
    if first_lib:
        e.to_csv(fname, sep='\t', index=False)
    else:
        e.to_csv(fname, sep='\t', index=False, header=None, mode='a')

    # library
    l = df.copy(deep=True)

    # create the alias
    l['aliases'] = 'ali-mortazavi:library_'+l.alias_assay+ \
                   '_'+l.alias_genotype+ \
                   '_'+l.tissue_desc+ \
                   '_'+l.age_desc+ \
                   '_'+l.model_organism_sex+ \
                   '_'+l.depth+ \
                   '_'+l.batch+ \
                   '_'+l.rep

    # # create the description
    # l['description'] = l.desc_assay+ \
    #                    ' '+l.desc_genotype+ \
    #                    ' '+l.tissue_desc+ \
    #                    ' '+l.age_desc+ \
    #                    ' '+l.model_organism_sex+ \
    #                    ' '+l.depth+ \
    #                    ' '+l.rep

    cols = ['aliases', 'biosample', 'nucleic_acid_term_name',
            'documents', 'construction_method', 'nucleic_acid_starting_quantity',
            'nucleic_acid_starting_quantity_units', 'strand_specificity',
            'average_fragment_size', 'lab', 'award', 'Sample ID']
    l = l[cols]
    l.drop_duplicates(inplace=True)

    # merge in alias with orig. df to have access for replicates, files
    l_alias = l[['aliases', 'Sample ID']]
    l_alias.rename({'aliases': 'library'}, axis=1, inplace=True)
    df = df.merge(l_alias, how='left', on='Sample ID')
    l.drop('Sample ID', axis=1, inplace=True)
    l.drop_duplicates(inplace=True)

    fname = 'library.tsv'
    if first_lib:
        l.to_csv(fname, sep='\t', index=False)
    else:
        l.to_csv(fname, sep='\t', index=False, header=None, mode='a')

    # replicate
    r = df.copy(deep=True)

    # create the alias
    r['aliases'] = 'ali-mortazavi:replicate_'+r.alias_assay+ \
                   '_'+r.alias_genotype+ \
                   '_'+r.tissue_desc+ \
                   '_'+r.age_desc+ \
                   '_'+r.model_organism_sex+ \
                   '_'+r.depth+ \
                   '_'+r.batch+ \
                   '_'+r.rep

    cols = ['aliases', 'library', 'experiment',
            'biological_replicate_number', 'technical_replicate_number',
            'Sample ID']
    r = r[cols]
    r.drop_duplicates(inplace=True)

    # merge in alias with orig. df to have access for files
    r_alias = r[['aliases', 'Sample ID']]
    r_alias.rename({'aliases': 'replicate'}, axis=1, inplace=True)
    df.drop('replicate', axis=1, inplace=True)
    df = df.merge(r_alias, how='left', on='Sample ID')
    r.drop('Sample ID', axis=1, inplace=True)
    r.drop_duplicates(inplace=True)

    fname = 'replicate.tsv'
    if first_lib:
        r.to_csv(fname, sep='\t', index=False)
    else:
        r.to_csv(fname, sep='\t', index=False, header=None, mode='a')

    if assay == 'illumina':

        # r1 fastqs
        r1 = df.copy(deep=True)

        # create the alias
        r1['aliases'] = 'ali-mortazavi:r1_fastq_'+r1.alias_assay+ \
                       '_'+r1.alias_genotype+ \
                       '_'+r1.tissue_desc+ \
                       '_'+r1.age_desc+ \
                       '_'+r1.model_organism_sex+ \
                       '_'+r1.depth+ \
                       '_'+r1.batch+ \
                       '_'+r1.rep

        # rename experiment and fname
        r1.rename({'experiment': 'dataset',
                   'fname': 'submitted_file_name'}, axis=1, inplace=True)

        # r1-specific things
        r1['file_format'] = 'fastq'
        r1['output_type'] = 'reads'
        r1['read_length'] = 115
        r1['run_type'] = 'single-ended'

        # limit only to r1
        r1 = r1.loc[~r1.submitted_file_name.str.contains('_R2.fastq.gz')]

        cols = ['aliases', 'submitted_file_name', 'dataset',
                'replicate', 'file_format', 'output_type',
                'platform', 'lab', 'award', 'read_length',
                'run_type', 'Sample ID']
        r1 = r1[cols]
        r1.drop_duplicates(inplace=True)

        # merge in alias with orig. df to have access for r1s
        r1_alias = r1[['aliases', 'Sample ID']]
        r1_alias.rename({'aliases': 'index_of'}, axis=1, inplace=True)
        df = df.merge(r1_alias, how='left', on='Sample ID')
        r1.drop('Sample ID', axis=1, inplace=True)
        r1.drop_duplicates(inplace=True)

        fname = 'r1_fastq.tsv'
        if first_lib:
            r1.to_csv(fname, sep='\t', index=False)
        else:
            r1.to_csv(fname, sep='\t', index=False, header=None, mode='a')

        # r2 fastqs
        r2 = df.copy(deep=True)

        # create the alias
        r2['aliases'] = 'ali-mortazavi:r2_fastq_'+r2.alias_assay+ \
                       '_'+r2.alias_genotype+ \
                       '_'+r2.tissue_desc+ \
                       '_'+r2.age_desc+ \
                       '_'+r2.model_organism_sex+ \
                       '_'+r2.depth+ \
                       '_'+r2.batch+ \
                       '_'+r2.rep

        # rename experiment and fname
        r2.rename({'experiment': 'dataset',
                   'fname': 'submitted_file_name'}, axis=1, inplace=True)

        # r2-specific things
        r2['file_format'] = 'fastq'
        r2['output_type'] = 'index reads'
        r2['read_length'] = 86
        r2['run_type'] = 'single-ended'
        r2['read_structure'] = '{"sequence_element": "UMI", "start": 1, "end": 10},{"sequence_element": "cell barcode", "start": 11, "end": 18},{"sequence_element": "cell barcode", "start": 49, "end": 56},{"sequence_element": "cell barcode", "start": 79, "end": 86}'

        # limit only to r2
        r2 = r2.loc[r2.submitted_file_name.str.contains('_R2.fastq.gz')]

        cols = ['aliases', 'dataset', 'submitted_file_name',
                'replicate', 'file_format', 'output_type',
                'platform', 'lab', 'award', 'read_length',
                'run_type', 'index_of', 'read_structure', 'Sample ID']
        r2 = r2[cols]

        # merge in alias with orig. df to have access for downstream files
        r2_alias = r2[['aliases', 'Sample ID']]
        r2_alias.rename({'aliases': 'r2_fastq'}, axis=1, inplace=True)
        df = df.merge(r2_alias, how='left', on='Sample ID')
        df.rename({'index_of': 'r1_fastq'}, axis=1, inplace=True)
        r2.drop('Sample ID', axis=1, inplace=True)
        r2.drop_duplicates(inplace=True)

        fname = 'r2_fastq.tsv'
        if first_lib:
            r2.to_csv(fname, sep='\t', index=False)
        else:
            r2.to_csv(fname, sep='\t', index=False, header=None, mode='a')

    first_lib = False

# finally just remove dupes from each thing
fname = 'biosample.tsv'
df = pd.read_csv(fname, sep='\t')
df.drop_duplicates(inplace=True)
df.to_csv(fname, sep='\t', index=False)

fname = 'experiment.tsv'
df = pd.read_csv(fname, sep='\t')
df.drop_duplicates(inplace=True)
df.to_csv(fname, sep='\t', index=False)

fname = 'library.tsv'
df = pd.read_csv(fname, sep='\t')
df.drop_duplicates(inplace=True)
df.to_csv(fname, sep='\t', index=False)

fname = 'replicate.tsv'
df = pd.read_csv(fname, sep='\t')
df.drop_duplicates(inplace=True)
df.to_csv(fname, sep='\t', index=False)

fname = 'r1_fastq.tsv'
df = pd.read_csv(fname, sep='\t')
df.drop_duplicates(inplace=True)
df.to_csv(fname, sep='\t', index=False)

fname = 'r2_fastq.tsv'
df = pd.read_csv(fname, sep='\t')
df.drop_duplicates(inplace=True)
df.to_csv(fname, sep='\t', index=False)
