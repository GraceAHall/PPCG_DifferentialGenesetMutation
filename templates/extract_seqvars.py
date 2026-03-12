
import re 
import pandas as pd 
import argparse
from copy import deepcopy
from extract_utils import CCFestimator

pd.set_option('display.max_rows', 500)
pd.set_option('display.max_columns', 500)
pd.set_option('display.width', 1000)
pd.set_option('display.max_colwidth', 80)
pd.options.display.float_format = '{:.2f}'.format


def main() -> None:
    args = load_cmdline_args()
    
    # load vcf
    df = load_vcf(args.vcf, args.vtype)
    print()
    print(df['func'].value_counts())
    
    # handle terms
    terms = select_terms(args)
    df = df[df['func'].isin(terms)].copy()
    print()
    print(df['func'].value_counts())
    
    # handle distances to genes
    if args.allow_intergenic_dist != 0:
        df = df[df['dist']<=args.allow_intergenic_dist].copy()
    # df = df.drop('dist', axis=1)
    print()
    print(df['func'].value_counts())

    # ccf estimation
    if df.shape[0] != 0:
        estimator = CCFestimator(purity=args.purity, scna_path=args.scna)
        df['est_ccf'] = df.apply(lambda row: estimator.est_ccf(row.chr, row.pos, row.VAF), axis=1)
    else:
        df = pd.DataFrame(columns=df.columns.to_list() + ['est_ccf'])

    # coordinates
    df['pos'] = df['pos'].astype(str)
    df['coords'] = df['chr'] + ':' + df['pos']
    df = df[['coords', 'REF', 'ALT', 'VAF', 'est_ccf', 'func', 'gene', 'dist']]
    print()
    print(df['func'].value_counts())
    
    df = df.rename(columns={'func': 'annotation'})
    df.to_csv(args.outfile, sep='\t', index=False, float_format='%.2f')

def load_vcf(filepath: str, vtype: str) -> pd.DataFrame:
    VCFFIELDS = ['chr','pos1','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT','NORMAL','TUMOUR']
    df = pd.read_csv(filepath, comment='#', sep='\t', header=None)
    df.columns = VCFFIELDS
    
    # no funny business
    df = df[df['FILTER']=='PASS'].copy()
    df['chr'] = df['chr'].astype(str)
    df['pos1'] = df['pos1'].astype(int)
    df = df[(df['chr'].str.isdigit()) | (df['chr'].isin(['X', 'Y']))].copy()
    
    # VAF
    assert vtype in ['SNV', 'INDEL']
    idx_tot_reads = 1 if vtype == 'INDEL' else 0
    idx_alt_reads = 2
    df['num_reads'] = df['TUMOUR'].apply(lambda x: x.split(':')[idx_tot_reads])
    df['var_reads'] = df['TUMOUR'].apply(lambda x: x.split(':')[idx_alt_reads])
    df['num_reads'] = df['num_reads'].astype(int)
    df['var_reads'] = df['var_reads'].astype(int)
    df['VAF'] = df['var_reads'] / df['num_reads']
    
    data = []
    for rec in df.itertuples():
        func_l, gene_l, dist_l = get_info(rec.INFO) # type: ignore
        data.append((
            rec.chr, rec.pos1, rec.REF, rec.ALT, rec.VAF,
            func_l, gene_l, dist_l
        ))

    df = pd.DataFrame.from_records(data, columns=['chr', 'pos', 'REF', 'ALT', 'VAF', 'func', 'gene', 'dist'])
    
    # bakuhatsu
    df = df.explode(['func', 'gene', 'dist']).copy()
    irows = df.shape[0]
    df = df.dropna().copy()
    frows = df.shape[0]
    print()
    print(f"dropped {irows-frows} rows with NA values for func, gene, or dist.")
    df['dist'] = df['dist'].astype(int)
    return df



ANNOVAR_TO_SO_LUT = {
    'UTR3': '3_prime_UTR_variant',
    'UTR5': '5_prime_UTR_variant',
    'intergenic': 'intergenic_variant',
    'intronic': 'intron_variant',
    'upstream': 'upstream_gene_variant',
    'downstream': 'downstream_gene_variant',
    'synonymous SNV': 'synonymous_variant',
    'unknown': 'unknown',
    'ncRNA_intronic': 'non_coding_transcript_intron_variant',
    'startloss': 'start_lost',
    'stopgain': 'stop_gained',
    'stoploss': 'stop_lost',
    'frameshift substitution': 'frameshift_variant',
    'nonsynonymous SNV': 'missense_variant',
    'splicing': 'splice_site_variant',
    'ncRNA_exonic': 'non_coding_transcript_exon_variant',
    'ncRNA_splicing': 'non_coding_transcript_splice_region_variant',
    'nonframeshift substitution': 'inframe_variant',
}
SO_SIMPLIFICATION_LUT = {
    'inframe_indel': 'inframe_variant',
    'inframe_deletion': 'inframe_variant',
    'inframe_insertion': 'inframe_variant',
    'splice_acceptor_variant': 'splice_site_variant',
    'splice_donor_variant': 'splice_site_variant',
    'splice_donor_5th_base_variant': 'splice_site_variant',
}

PATTERN_FUNC = r'.*Func\.refGene=(.+?(?=;Gene)).*'
PATTERN_EXFUNC = r'.*ExonicFunc\.refGene=(.+?(?=;AAChange)).*'
PATTERN_GENE = r'.*Gene\.refGene=(.+(?=;GeneDetail)).*'
PATTERN_DIST = r'.*GeneDetail.refGene=((dist=\w+;)+).*'

def get_info(info: str):
    func = re.match(PATTERN_FUNC, info).group(1) # type: ignore
    func_l = func.split(';')
    func_l = [ANNOVAR_TO_SO_LUT.get(x, x) for x in func_l]
    func_l = [SO_SIMPLIFICATION_LUT.get(x, x) for x in func_l]

    m = re.match(PATTERN_EXFUNC, info)
    if m is not None:
        exfunc_l = [m.group(1)]
    else:
        exfunc_l = []
    exfunc_l = [ANNOVAR_TO_SO_LUT.get(x, x) for x in exfunc_l]
    exfunc_l = [SO_SIMPLIFICATION_LUT.get(x, x) for x in exfunc_l]

    gene = re.match(PATTERN_GENE, info).group(1) # type: ignore
    gene_l = gene.split(';')

    m = re.match(PATTERN_DIST, info)
    if m is not None:
        dist_l = m.group(1).rstrip(';').split(';')
        dist_l = [x.replace('dist=', '') for x in dist_l]
    else:
        dist_l = []

    # swap out funcs for exonic funcs
    if len(exfunc_l) > 0:
        # assert len(func_l) == 1
        # assert len(set(gene_l)) == 1
        assert len(exfunc_l) == 1
        assert len(dist_l) == 0
        func_l_new = []
        ptr = 0
        for func in func_l:
            if func == 'exonic':
                func_l_new.append(exfunc_l[ptr])
                if ptr != len(exfunc_l)-1:
                    ptr += 1 
            else:
                func_l_new.append(func)
        func_l = func_l_new

    if len(func_l) != len(gene_l) or len(gene_l) != len(dist_l):
        func_l, gene_l, dist_l = handle_annotation_weirdness(func_l, gene_l, dist_l)

    out_funcs = []
    out_genes = []
    out_dists = []
    for func, gene, dist in zip(func_l, gene_l, dist_l):
        if gene == 'NONE' and dist == 'NONE':
            continue
        if gene == 'NONE' or dist == 'NONE':
            print()
            print('[WARN]: unresolvable ANNOVAR annotation.')
            print(f"funcs={func_l}")
            print(f"genes={gene_l}")
            print(f"dists={dist_l}")
            raise NotImplementedError
        out_funcs.append(func)
        out_genes.append(gene)
        out_dists.append(dist)

    return out_funcs, out_genes, out_dists

def handle_annotation_weirdness(func_l: list[str], gene_l: list[str], dist_l: list[str]):
    intergenic = ['intergenic_variant', 'downstream_gene_variant', 'downstream_gene_variant']
    # print()
    # print(func_l, gene_l, dist_l)

    if len(func_l) != len(gene_l):
        # attempt to resolve situation where a 'NONE' has been injected for some reason.
        # god I hate ANNOVAR.
        gene_l_nn = [x for x in gene_l if x!='None']
        if len(func_l) == len(gene_l_nn):
            gene_l = gene_l_nn
        # attempt to resolve situation where a single function seems to correspond to all genes. 
        elif len(func_l)==1 and len(gene_l)>1:
            func_l = func_l * len(gene_l)
        elif len(gene_l) > len(func_l):
            print()
            print('[WARN]: unresolvable ANNOVAR annotation.')
            print(f"funcs={func_l}")
            print(f"genes={gene_l}")
            print(f"dists={dist_l}")
            while len(func_l) < len(gene_l):
                func_l.append('unknown')
        else:
            raise NotImplementedError

    if len(func_l) != len(gene_l):
        raise NotImplementedError
    
    ptr = 0
    dist_l_new = []
    for func, gene in zip(func_l, gene_l):
        if func in intergenic:
            dist_l_new.append(dist_l[ptr])
            if ptr != len(dist_l)-1:
                ptr += 1 
        else:
            dist_l_new.append(0)
    # make sure all distances have been consumed
    if len(dist_l) > 0:
        if ptr!=len(dist_l)-1:
            print(ptr)
            print(f"funcs={func_l}")
            print(f"genes={gene_l}")
            print(f"dists={dist_l}")
            raise NotImplementedError
        assert ptr==len(dist_l)-1
    dist_l = dist_l_new
    
    if len(func_l) != len(dist_l):
        raise NotImplementedError
    
    return func_l, gene_l, dist_l
    

SO_TERMS = [
    'missense_variant',
    'inframe_variant',
    'frameshift_variant',
    'start_lost',
    'stop_gained',
    'stop_lost',
    'splice_site_variant',
    '3_prime_UTR_variant',
    '5_prime_UTR_variant',

    'non_coding_transcript_exon_variant',
    'non_coding_transcript_splice_region_variant',

    'upstream_gene_variant',
    'downstream_gene_variant',
    'intergenic_variant',
    # 'synonymous_variant',
    # 'intron_variant',
    # 'non_coding_transcript_intron_variant',
    'unknown',
]

def select_terms(args: argparse.Namespace) -> list[str]:
    terms = deepcopy(SO_TERMS)
    if not args.allow_utr:
        terms = [t for t in terms if 'UTR' not in t]
    if not args.allow_splice:
        terms = [t for t in terms if 'splice' not in t]
    if not args.allow_noncoding:
        terms = [t for t in terms if 'non_coding' not in t]
    if args.allow_intergenic_dist == 0:
        banned = ['intergenic_variant', 'upstream_gene_variant', 'downstream_gene_variant']
        terms = [t for t in terms if t not in banned]
    # maybe TODO? 
    # if args.allow_all:
    #     return terms
    # if not args.allow_intron:
    #     terms = [t for t in terms if 'intron_variant' not in t]
    # # if not args.allow_synonymous:
    # #     terms = [t for t in terms if t!='synonymous_variant']
    return terms 


def load_cmdline_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description='Extracts SNVs/INDELs from a PPCG donor VCF file.')
    parser.add_argument('--vcf', type=str, required=True, help='Path to VCF file. ')
    parser.add_argument('--scna', type=str, required=True, help='Path to Battenberg SCNA segments (for CCF estimation).')
    parser.add_argument('--purity', type=float, required=True, help='Sample purity/cellularity (for CCF estimation).')
    parser.add_argument('--vtype', required=True, choices=['SNV', 'INDEL'], help='Variant type')
    parser.add_argument('--allow-utr', action='store_true', help='Include 3\' and 5\' UTR variants.')
    parser.add_argument('--allow-splice', action='store_true', help='Include splice variants.')
    parser.add_argument('--allow-noncoding', action='store_true', help='Include non-coding variants.')
    parser.add_argument('--allow-intergenic-dist', type=int, default=0, help='Include upsteam/downstream variants within <dist> of gene.')
    parser.add_argument('--outfile', type=str, required=True, help='Path to output file containing extracted variants')
    # parser.add_argument('--allow-all', action='store_true', help='Include all variants.')
    # parser.add_argument('--allow-synonymous', action='store_true', help='Include synonymous exon variants.')
    # parser.add_argument('--allow-intron', action='store_true', help='Include intron variants.')
    args = parser.parse_args()
    return args

if __name__ == '__main__':
    main()




