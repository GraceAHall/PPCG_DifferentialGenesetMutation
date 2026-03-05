
import re 
import sys 
import argparse
import warnings
import numpy as np
import pandas as pd 
import matplotlib.pyplot as plt
import seaborn as sns

warnings.filterwarnings('ignore')

pd.options.display.float_format = "{:,.3f}".format
pd.set_option('display.max_colwidth', None)
pd.set_option('display.max_columns', None)
pd.set_option('display.width', 200)


def main() -> None:
    args = load_cmdline_args()
    id2type, id2curr, curr2id, alt2id = load_hgnc(args)
    
    # load & standardise mutations file
    muts = pd.read_csv(args.mutations, sep='\t', header=0)
    muts = muts[muts['annotation']!='CNApeak'].copy()
    muts = standardise_hgnc_muts(muts, id2type=id2type, id2curr=id2curr, curr2id=curr2id, alt2id=alt2id)

    # load & standardise genesets file
    gframe = pd.read_csv(args.genesets, sep='\t', header=0)
    gframe = standardise_hgnc_genesets(gframe, id2type=id2type, id2curr=id2curr, curr2id=curr2id, alt2id=alt2id)
    
    # load gff, est gene sizes, standardise size information symbols
    valid = select_gff_transcripts(args)
    gff = load_gff_gencode(args, valid)
    sframe = generate_sizes_table(gff)
    sframe = standardise_hgnc_sizes(sframe, id2type=id2type, id2curr=id2curr, curr2id=curr2id, alt2id=alt2id)

    # filter non-protein-coding genes
    sframe, gframe, muts = filter_non_protein_coding_genes(sframe, gframe, muts)

    # ensure missing mutated / geneset genes have a size. 
    sframe = augment_sizes_table(sframe, gframe, muts)

    # write to disk
    write_files(args, muts, gframe, sframe)



##################################################################################
### Load GENCODE v31lift37 gene annotations -> Calculate gene size information ###
##################################################################################
# Juri states that ANNOVAR 2020Jun08 edition was used for SNV/INDEL annotations.
# In theory, ANNOVAR should be using this version of GENCODE for it's annotations. 
# Therefore all SNV/INDEL mutation symbols should match this GENCODE GFF file.

def _get_gencode_name(text: str) -> str:
    PATTERN = r'^.*?gene_name=([A-Za-z0-9_\.-]+).*?$'
    m = re.match(PATTERN, text)
    assert m is not None
    return m.group(1)

def _get_gencode_gtype(text: str) -> str:
    PATTERN = r'^.*?gene_type=([A-Za-z0-9_-]+).*?$'
    m = re.match(PATTERN, text)
    assert m is not None
    return m.group(1)

def _get_gencode_ttype(text: str) -> str:
    PATTERN = r'^.*?transcript_type=([A-Za-z0-9_-]+).*?$'
    m = re.match(PATTERN, text)
    assert m is not None
    return m.group(1)

def _get_gencode_tid(text: str) -> str:
    PATTERN = r'^.*;transcript_id=([A-Z0-9\._]+).*$'
    m = re.match(PATTERN, text)
    assert m is not None
    return m.group(1)

def _get_gencode_tags(text: str) -> str:
    PATTERN = r'^.*;tag=([^;]+).*$'
    m = re.match(PATTERN, text)
    assert m is not None
    return m.group(1)

def _get_gencode_hgnc_id(text: str) -> None|str:
    PATTERN = r'^.*?hgnc_id=HGNC:([0-9]+).*?$'
    m = re.match(PATTERN, text)
    if m is None:
        return None
    return m.group(1)

def select_gff_transcripts(args: argparse.Namespace) -> set[str]:
    
    banned_tags = [
        'non_canonical_TEC',
        'non_canonical_U12',
        'non_canonical_conserved',
        'non_canonical_genome_sequence_error',
        'non_canonical_other',
        'non_canonical_polymorphism',
        'non_submitted_evidence',
        'not_organism_supported',
        'not_best_in_genome_evidence',
        'readthrough_transcript',
        'RNA_Seq_supported_only',
        'retained_intron_CDS',
        'retained_intron_final',
        'retained_intron_first',
    ]
    
    # # additionally remove these gff genes cuz they're bent
    # banned_transcripts = [
    #     'ENSG00000033050.9_4',
    #     'ENSG00000286169.1_1',
    #     'ENSG00000163635.18_4',
    #     'ENSG00000145075.13_7',
    #     'ENSG00000261480.1_6',
    #     'ENSG00000184047.19_6',
    #     'ENSG00000145194.18_3',
    #     'ENSG00000286070.1_1',
    #     'ENSG00000284024.2_4',
    #     'ENSG00000254093.9_6',
    #     'ENSG00000168255.20_5',
    #     'ENSG00000285053.1_3',
    # ]
    # sframe = sframe[~sframe['gene_id'].isin(banned_transcripts)].copy()

    out = set()
    with open(args.gff, 'r') as fp:

        line = fp.readline().strip('\n')
        while line:
            if not line.startswith('#'):
                lsplit = line.split('\t')
                if lsplit[2] == 'transcript':
                    ttype = _get_gencode_ttype(lsplit[8])
                    if ttype == 'protein_coding':
                        tid = _get_gencode_tid(lsplit[8])
                        tags = _get_gencode_tags(lsplit[8])
                        if not any([x in tags for x in banned_tags]):
                            out.add(tid)
            line = fp.readline().strip('\n')
    
    return out 

def load_gff_gencode(args: argparse.Namespace, valid: set) -> pd.DataFrame:
    
    data = []
    with open(args.gff, 'r') as fp:
        i = 0
        
        line = fp.readline().strip('\n')
        while line:
            if line.startswith('#'):
                line = fp.readline().strip('\n')
                continue 
            
            i += 1
            if i % 10_000 == 0:
                print(f"processed {i} records...", end='\r')
            
            lsplit = line.split('\t')
            ftype = str(lsplit[2])

            if ftype!='exon' and ftype!='CDS':
                line = fp.readline().strip('\n')
                continue

            tid = _get_gencode_tid(lsplit[8])
            if tid not in valid:
                line = fp.readline().strip('\n')
                continue

            gene = _get_gencode_name(lsplit[8])
            hgnc_id = _get_gencode_hgnc_id(lsplit[8])
            chrom = lsplit[0].replace('chr', '')
            start = int(lsplit[3])
            end = int(lsplit[4])
            data.append((gene, hgnc_id, ftype, chrom, start, end))

            line = fp.readline().strip('\n')

    print(f"processed {i} records... done")
    df = pd.DataFrame.from_records(data, columns=['gene', 'hgnc_id', 'feature', 'chrom', 'start', 'end'])
    df['gene'] = df['gene'].astype(str)
    # df['tid'] = df['tid'].astype(str)
    df['feature'] = df['feature'].astype(str)
    df['chrom'] = df['chrom'].astype(str)
    df['start'] = df['start'].astype(int)
    df['end'] = df['end'].astype(int)
    df['span'] = df['end'] - df['start']
    df['span'] = df['span'].astype(int)
    df = df.sort_values('gene')
    return df 


############################
### Calculate Gene Sizes ###
############################

def generate_sizes_table(gff: pd.DataFrame):
    n_genes = gff['gene'].nunique()
    hid_LUT = gff[gff['hgnc_id'].notna()].drop_duplicates('gene').set_index('gene')['hgnc_id'].to_dict()
    data = []
    i = 0

    grouped_chunks = {name: group for name, group in gff.groupby('gene')}
    for gene, df in grouped_chunks.items():
        i += 1
        if i % 100 == 0:
            print(f"processed {i}/{n_genes} genes...", end='\r')
        echunk = df[df['feature']=='exon']
        cchunk = df[df['feature']=='CDS']
        if echunk.shape[0] == 0:
            print(f"Skipped {gene} (no exon information).")
            continue
        if cchunk.shape[0] == 0:
            print(f"Skipped {gene} (no CDS information).")
            continue
        
        data.append((
            gene, 
            hid_LUT.get(gene, None),
            _get_multichrom_gene_span(echunk),
            _get_logical_or_segments_len(echunk),
            _get_logical_or_segments_len(cchunk)
        ))
    print(f"processed {i}/{n_genes} genes... done.")

    df = pd.DataFrame.from_records(data, columns=['gene', 'hgnc_id', 'span', 'cum_exon_len', 'cum_cds_len'])
    df['span'] = df['span'].astype(int)
    df['cum_exon_len'] = df['cum_exon_len'].astype(int)
    df['cum_cds_len'] = df['cum_cds_len'].astype(int)
    return df

def _get_logical_or_segments_len(df: pd.DataFrame) -> int:
    chrom_sums = []
    chrom_dfs = []
    if df['chrom'].nunique() > 1:
        for chrom in df['chrom'].unique():
            chrom_dfs.append(df[df['chrom']==chrom].copy())
    elif df['chrom'].nunique() == 1:
        chrom_dfs.append(df)
    else:
        raise NotImplementedError
    
    for df in chrom_dfs:
        left = df['start'].min()
        right = df['end'].max()
        offset = left
        width = right - left
        mask = np.zeros(width+1)
        for rec in df.itertuples():
            lb = rec.start - offset
            ub = rec.end - offset
            mask[lb:ub] = 1
        chrom_sums.append(mask.sum())
            
    return sum(chrom_sums)

def _get_multichrom_gene_span(df: pd.DataFrame) -> int:
    cumspan = 0
    chrom_dfs = []

    if df['chrom'].nunique() > 1:
        for chrom in df['chrom'].unique():
            chrom_dfs.append(df[df['chrom']==chrom].copy())
    elif df['chrom'].nunique() == 1:
        chrom_dfs.append(df)
    else:
        raise NotImplementedError
    
    for df in chrom_dfs:
        left = df['start'].min()
        right = df['end'].max()
        cumspan += right - left 
    return cumspan



######################
### Load HGNC LUTs ###
######################


def load_hgnc(args: argparse.Namespace):
    collisions = 0
    id2type = dict()
    id2curr = dict()
    curr2id = dict()
    alt2id = dict()

    with open(args.hgnc, 'r') as fp:
        line = fp.readline().strip('\n')
        line = fp.readline().strip('\n')
        while line:
            lsplit = line.split('\t')
            locus_group = lsplit[3]
            curr_id = lsplit[0].replace('HGNC:', '')
            curr_sym = lsplit[1]
            alias_syms = lsplit[8].strip('"').split('|')
            prev_syms = lsplit[10].strip('"').split('|')
            id2type[curr_id] = locus_group
            id2curr[curr_id] = curr_sym
            curr2id[curr_sym] = curr_id
            for alt_sym in alias_syms + prev_syms:
                if alt_sym == '':
                    continue 
                if alt_sym in alt2id and alt2id[alt_sym]!=curr_id:
                    collisions += 1
                alt2id[alt_sym] = curr_id
            line = fp.readline().strip('\n')
    print(len(curr2id))
    print(len(alt2id))
    print(collisions)
    print()
    return id2type, id2curr, curr2id, alt2id



###################################
### Standardise gene sizes data ###
###################################

def fetch_hgnc_id(text: str, alt2id: dict, curr2id: dict) -> str|None:
    if text in curr2id:
        return curr2id[text]
    elif text in alt2id:
        return alt2id[text]
    return None

def standardise_hgnc_sizes(table: pd.DataFrame, id2type: dict, id2curr: dict, curr2id: dict, alt2id: dict) -> pd.DataFrame:
    df = table.copy()
    
    # ids
    mask = df['hgnc_id'].isna()
    df.loc[mask, 'ID'] = df.loc[mask, 'gene'].apply(fetch_hgnc_id, alt2id=alt2id, curr2id=curr2id)
    df.loc[~mask, 'ID'] = df.loc[~mask, 'hgnc_id']

    # syms 
    mask = df['ID'].notna()
    df.loc[mask, 'SYM'] = df.loc[mask, 'ID'].apply(lambda x: id2curr.get(x, None))
    
    # functions
    mask = df['ID'].notna()
    df.loc[mask, 'FUNC'] = df.loc[mask, 'ID'].apply(lambda x: id2type.get(x, None))

    print(f"{df.shape} records")
    print(f"{df['gene'].notna().sum()} with 'gene'")
    print(f"{df['hgnc_id'].notna().sum()} with 'hgnc_id'")
    print(f"{df['ID'].notna().sum()} with 'ID'")
    print(f"{df['SYM'].notna().sum()} with 'SYM'")
    print(f"{df['FUNC'].notna().sum()} with 'FUNC'")
    print(f"{(df['ID']!=df['hgnc_id']).sum()} mismatch ids")
    print(f"{(df['gene']!=df['SYM']).sum()} mismatch gene symbols")
    print()
    print(df['FUNC'].value_counts(dropna=False))

    return df 


########################################
### Load & Standardise mutation data ###
########################################

def standardise_hgnc_muts(table: pd.DataFrame, id2type: dict, id2curr: dict, curr2id: dict, alt2id: dict) -> pd.DataFrame:
    df = table.copy()
    
    # ids
    df['ID'] = df['gene'].apply(fetch_hgnc_id, alt2id=alt2id, curr2id=curr2id)

    # syms & funcs
    mask = df['ID'].notna()
    df.loc[mask, 'SYM'] = df.loc[mask, 'ID'].apply(lambda x: id2curr.get(x, None))
    df.loc[mask, 'FUNC'] = df.loc[mask, 'ID'].apply(lambda x: id2type.get(x, None))

    temp = df.drop_duplicates('gene')
    print(f"{temp.shape} records")
    print(f"{temp['ID'].notna().sum()} with 'ID'")
    print(f"{temp['SYM'].notna().sum()} with 'SYM'")
    print(f"{temp['FUNC'].notna().sum()} with 'FUNC'")
    print(f"{(temp['gene']!=temp['SYM']).sum()} mismatch gene symbols")
    print()
    print('variants missing HGNC ID:')
    print(temp[temp['ID'].isna()]['vclass'].value_counts())

    return df 


#######################################
### Load & Standardise geneset data ###
#######################################

def standardise_hgnc_genesets(table: pd.DataFrame, id2type: dict, id2curr: dict, curr2id: dict, alt2id: dict) -> pd.DataFrame:
    df = table.copy()
    
    # ids
    df['ID'] = df['Gene'].apply(fetch_hgnc_id, alt2id=alt2id, curr2id=curr2id)
    
    # syms & funcs
    mask = df['ID'].notna()
    df.loc[mask, 'SYM'] = df.loc[mask, 'ID'].apply(lambda x: id2curr.get(x, None))
    df.loc[mask, 'FUNC'] = df.loc[mask, 'ID'].apply(lambda x: id2type.get(x, None))

    temp = df.drop_duplicates('Gene')
    print(f"{temp.shape} records")
    print(f"{temp['ID'].notna().sum()} with 'ID'")
    print(f"{temp['SYM'].notna().sum()} with 'SYM'")
    print(f"{temp['FUNC'].notna().sum()} with 'FUNC'")
    print(f"{(temp['Gene']!=temp['SYM']).sum()} mismatch gene symbols")
    print()

    return df 


#######################################
### Filter non-protein-coding genes ###
#######################################

def filter_non_protein_coding_genes(sframe: pd.DataFrame, gframe: pd.DataFrame, muts: pd.DataFrame):
    for label, df in zip(['sizes', 'genesets', 'mutations'], [sframe, gframe, muts]):
        print()
        print(f"--- {label} ---")
        print(f"{df.shape[0]} total records.")
        print(f"{len(df['ID'].unique())} known genes.")
        # print(df.groupby('FUNC')['ID'].nunique())
        df = df[df['FUNC']=='protein-coding gene']
        print(f"{len(df['ID'].unique())} protein-coding genes will be kept.")

    muts = muts[muts['FUNC']=='protein-coding gene'].copy()
    sframe = sframe[sframe['FUNC']=='protein-coding gene'].copy()
    gframe = gframe[gframe['FUNC']=='protein-coding gene'].copy()

    assert muts['ID'].isna().sum() == 0
    assert sframe['ID'].isna().sum() == 0
    assert gframe['ID'].isna().sum() == 0

    return sframe, gframe, muts



##################################################################
### check overlap & add missing genes to size data (as median) ###
##################################################################

def augment_sizes_table(sframe: pd.DataFrame, gframe: pd.DataFrame, muts: pd.DataFrame):
    # first, drop non-required cols from sframe
    sframe = sframe[['SYM', 'span', 'cum_exon_len', 'cum_cds_len']].copy()

    # calc imputed values
    med_span = sframe['span'].median()
    med_exon = sframe['cum_exon_len'].median()
    med_cds = sframe['cum_cds_len'].median()
    
    # identify missing genes
    size_genes = set(sframe['SYM'].unique())
    gset_genes = set(gframe['SYM'].unique())
    mut_genes = set(muts['SYM'].unique())
    targets = set()
    targets = targets | (mut_genes-size_genes)
    targets = targets | (gset_genes-size_genes)

    print()
    print(f"size genes: {len(size_genes)}")
    print(f"geneset genes: {len(gset_genes)}")
    print(f"mutated genes: {len(mut_genes)}")
    print(f"geneset genes without size info: {len(gset_genes-size_genes)}")
    print(f"mutated genes without size info: {len(mut_genes-size_genes)}")
    print(f"total genes to add: {len(targets)}")

    data = []
    for sym in sorted(targets):
        data.append((sym, med_span, med_exon, med_cds))
    if len(data) > 0:
        df = pd.DataFrame.from_records(data, columns=sframe.columns.to_list())
        sframe = pd.concat([sframe, df], ignore_index=True)
    sframe['span'] = sframe['span'].astype(int)
    sframe['cum_exon_len'] = sframe['cum_exon_len'].astype(int)
    sframe['cum_cds_len'] = sframe['cum_cds_len'].astype(int)
    return sframe

#####################
### write to file ###
#####################

def write_files(args: argparse.Namespace, muts: pd.DataFrame, gframe: pd.DataFrame, sframe: pd.DataFrame):

    # mutations
    df = muts.drop(['gene', 'ID', 'FUNC'], axis=1)
    df = df.rename(columns={'SYM': 'gene'})
    df.to_csv(args.outfile_muts, sep='\t', index=False, float_format='%.3f')

    # genesets
    df = gframe[['SYM', 'Pathway']].copy()
    df = df.rename(columns={'SYM': 'gene', 'Pathway': 'geneset'})
    df.to_csv(args.outfile_gsets, sep='\t', index=False)

    # sizes
    df = sframe.copy()
    df = df.rename(columns={'SYM': 'gene'})
    df.to_csv(args.outfile_sizes, sep='\t', index=False)



def load_cmdline_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description='Standardises symbols.')
    parser.add_argument('--mutations', type=str, required=True, help='Path to extracted mutations file (tsv).')
    parser.add_argument('--genesets', type=str, required=True, help='Path to GSEA genesets file (tsv).')
    parser.add_argument('--gff', type=str, required=True, help='Path to GENCODE v31lift37.basic annotations file.')
    parser.add_argument('--hgnc', type=str, required=True, help='Path to HGNC file (txt).')
    parser.add_argument('--outfile-muts', type=str, required=True, help='Path to processed mutations outfile.')
    parser.add_argument('--outfile-gsets', type=str, required=True, help='Path to processed genesets outfile.')
    parser.add_argument('--outfile-sizes', type=str, required=True, help='Path to generated sizes outfile.')

    args = parser.parse_args()
    return args

if __name__ == '__main__':
    main()



    
# def load_gff_gencode(args: argparse.Namespace) -> pd.DataFrame:

#     def get_name(text: str) -> str:
#         PATTERN = r'^.*?gene_name=([A-Za-z0-9_\.-]+).*?$'
#         m = re.match(PATTERN, text)
#         assert m is not None
#         return m.group(1)

#     def get_gtype(text: str) -> str:
#         PATTERN = r'^.*?gene_type=([A-Za-z0-9_-]+).*?$'
#         m = re.match(PATTERN, text)
#         assert m is not None
#         return m.group(1)
    
#     def get_tid(text: str) -> str:
#         PATTERN = r'^.*;transcript_id=([A-Z0-9\._]+).*$'
#         m = re.match(PATTERN, text)
#         assert m is not None
#         return m.group(1)
    
#     def get_gid(text: str) -> str:
#         PATTERN = r'^.*;gene_id=([A-Z0-9\._-]+).*$'
#         m = re.match(PATTERN, text)
#         assert m is not None
#         return m.group(1)
    
#     def get_ttype(text: str) -> str:
#         PATTERN = r'^.*?transcript_type=([A-Za-z0-9_-]+).*?$'
#         m = re.match(PATTERN, text)
#         assert m is not None
#         return m.group(1)
    
#     def get_tags(text: str) -> str|None:
#         PATTERN = r'^.*;tag=([^;]+).*$'
#         m = re.match(PATTERN, text)
#         if m is not None:
#             return m.group(1)
#         return None

#     data = []
#     with open(args.gff, 'r') as fp:
#         line = fp.readline().strip('\n')
#         i = 0
#         while line:
#             if line.startswith('#'):
#                 line = fp.readline().strip('\n')
#                 continue 
            
#             i += 1
#             if i % 10_000 == 0:
#                 print(f"processed {i} records...", end='\r')
            
#             lsplit = line.split('\t')
#             chrom = lsplit[0].replace('chr', '')

#             if lsplit[2] != 'transcript':
#                 line = fp.readline().strip('\n')
#                 continue
            
#             start = int(lsplit[3])
#             end = int(lsplit[4])
#             gene = get_name(lsplit[8])
#             gtype = get_gtype(lsplit[8])
#             ttype = get_ttype(lsplit[8])
#             gid = get_gid(lsplit[8])
#             tid = get_tid(lsplit[8])
#             tags = get_tags(lsplit[8])

#             if gtype == 'protein_coding':
#                 data.append((gene, gid, tid, ttype, chrom, start, end, end-start, tags))

#             line = fp.readline().strip('\n')

#     print(f"processed {i} records... done")
#     df = pd.DataFrame.from_records(data, columns=['gene', 'gid', 'tid', 'ttype', 'chrom', 'start', 'end', 'span', 'tags'])
#     df['chrom'] = df['chrom'].astype(str)
#     df['start'] = df['start'].astype(int)
#     df['end'] = df['end'].astype(int)
#     df['span'] = df['span'].astype(int)
#     df['gene'] = df['gene'].astype(str)

#     banned = [
#         'non_canonical_TEC',
#         'non_canonical_U12',
#         'non_canonical_conserved',
#         'non_canonical_genome_sequence_error',
#         'non_canonical_other',
#         'non_canonical_polymorphism',
#         'non_submitted_evidence',
#         'not_organism_supported',
#         'not_best_in_genome_evidence',
#         'readthrough_transcript',
#         'RNA_Seq_supported_only',
#         'retained_intron_CDS',
#         'retained_intron_final',
#         'retained_intron_first',
#     ]

#     df['ignore'] = df['tags'].apply(lambda tag: any([x in tag for x in banned]))
#     df = df[df['ignore']==False].drop('ignore', axis=1).copy()
#     return df


# def get_logical_or_segments_len(df: pd.DataFrame) -> int:
#     chrom_sums = []
#     chrom_dfs = []
#     if df['chrom'].nunique() > 1:
#         for chrom in df['chrom'].unique():
#             chrom_dfs.append(df[df['chrom']==chrom].copy())
#     elif df['chrom'].nunique() == 1:
#         chrom_dfs.append(df)
#     else:
#         raise NotImplementedError
    
#     for df in chrom_dfs:
#         left = df['start'].min()
#         right = df['stop'].max()
#         offset = left
#         width = right - left
#         mask = np.zeros(width+1)
#         for rec in df.itertuples():
#             lb = rec.start - offset
#             ub = rec.stop - offset
#             mask[lb:ub] = 1
#         chrom_sums.append(mask.sum())
            
#     return sum(chrom_sums)

# def generate_sizes_table(gff: pd.DataFrame):
#     all_genes = sorted(list(gff['SYM'].unique()))
#     sframe = pd.DataFrame(index=all_genes)

#     # span
#     temp = gff[gff['feature']=='gene'].sort_values('span', ascending=False).drop_duplicates('SYM')
#     sframe['span'] = temp.set_index('SYM')['span']

#     # cum_exon / cum_cds
#     i = 0
#     grouped_chunks = {name: group for name, group in gff.groupby('SYM')}
#     for gene, df in grouped_chunks.items():
#         i += 1
#         if i % 100 == 0:
#             print(f"processed {i}/{len(all_genes)} genes...", end='\r')
#             sframe = sframe.copy()
#         echunk = df[df['feature']=='exon']
#         cchunk = df[df['feature']=='CDS']
#         if echunk.shape[0] == 0:
#             print(f"Skipped {gene} (no exon information).")
#             continue 
#         if cchunk.shape[0] == 0:
#             print(f"Skipped {gene} (no CDS information).")
#             continue 
#         sframe.loc[gene, 'cum_exon_len'] = get_logical_or_segments_len(echunk)
#         sframe.loc[gene, 'cum_cds_len'] = get_logical_or_segments_len(cchunk)
#     print(f"processed {i}/{len(all_genes)} genes... done.")

#     mask = (sframe['cum_exon_len'].isna()) | (sframe['cum_cds_len'].isna())
#     sframe = sframe.loc[~mask].copy()
#     sframe = sframe.astype(int)
#     sframe = sframe.reset_index()
#     sframe = sframe.rename(columns={'index': 'SYM'})
#     return sframe


# def standardise_hgnc_gff(table: pd.DataFrame, id2type: dict, id2curr: dict, curr2id: dict, alt2id: dict) -> pd.DataFrame:
#     df = table.copy()
    
#     # ids
#     mask = df['hgnc_id'].isna()
#     df.loc[mask, 'ID'] = df.loc[mask, 'gene'].apply(fetch_hgnc_id, alt2id=alt2id, curr2id=curr2id)
#     df.loc[~mask, 'ID'] = df.loc[~mask, 'hgnc_id']

#     # syms 
#     mask = df['ID'].notna()
#     df.loc[mask, 'SYM'] = df.loc[mask, 'ID'].apply(lambda x: id2curr.get(x, None))
    
#     # functions
#     mask = df['ID'].notna()
#     df.loc[mask, 'FUNC'] = df.loc[mask, 'ID'].apply(lambda x: id2type.get(x, None))

#     temp = df.drop_duplicates('gene')
#     print(f"{temp.shape} records")
#     print(f"{temp['gene'].notna().sum()} with 'gene'")
#     print(f"{temp['hgnc_id'].notna().sum()} with 'hgnc_id'")
#     print(f"{temp['ID'].notna().sum()} with 'ID'")
#     print(f"{temp['SYM'].notna().sum()} with 'SYM'")
#     print(f"{temp['FUNC'].notna().sum()} with 'FUNC'")
#     print(f"{(temp['ID']!=temp['hgnc_id']).sum()} mismatch ids")
#     print(f"{(temp['gene']!=temp['SYM']).sum()} mismatch gene symbols")

#     return df 