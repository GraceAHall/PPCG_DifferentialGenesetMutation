
import re 
import sys
import argparse
import pandas as pd 
from abc import ABC, abstractmethod

pd.set_option('display.max_rows', 500)
pd.set_option('display.max_columns', 500)
pd.set_option('display.width', 1000)
pd.set_option('display.max_colwidth', 40)
pd.options.display.float_format = '{:.2f}'.format

from extract_utils import load_scna

def main() -> None:
    args = load_cmdline_args()
    segs = load_scna(args.scna, allow_subclonal=args.allow_subclonal)
    
    # print()
    # print(segs.shape)
    # print(segs.head(10))
    
    # gff = load_gff_refseq(args.gff, args.allow_noncoding, args.allow_x)
    gff = load_gff_gencode(args.gff, args.allow_noncoding, args.allow_x)

    # print()
    # print(gff.shape)
    # print(gff.head(10))
    
    # table = extract_refseq(segs, gff, args)
    table = extract_gencode(segs, gff, args)
    table.to_csv(args.outfile, sep='\t', index=False)


###############
### FILE IO ###
###############

def load_gff_refseq(filepath: str, allow_noncoding: bool, allow_x: bool) -> pd.DataFrame:

    def get_chrom(seqid: str) -> str:
        lut = {23: 'X', 24: 'Y'}
        suffix = int(seqid.split('.')[0][-2:])
        return lut.get(suffix, str(suffix))
    
    def get_gene(info: str) -> str:
        pattern = r'.*gene=([^;]+);.*'
        m = re.match(pattern, info)
        if m is None:
            print(info)
            raise NotImplementedError
        return m.group(1)
    
    df = pd.read_csv(filepath, sep='\t', comment='#')
    df.columns = ['seqid','source','type','start','end','score','strand','phase','info']
    
    # remove annotations we're not interested in
    banned = [
        'region', 
        'match', 'cDNA_match', 
        'C_gene_segment', 'V_gene_segment', 
        'telomerase_RNA', 'RNase_MRP_RNA', 'RNase_P_RNA', 
        'transcript', 'exon', 'CDS', 'mRNA'
    ]
    df = df[~df['type'].isin(banned)].copy()

    # remove patches, Y and X
    df = df[df['seqid'].str.startswith('NC_')].copy()
    df['chr'] = df['seqid'].apply(get_chrom)
    df = df[df['chr']!='Y'].copy()
    if not allow_x:
        df = df[df['chr']!='X'].copy()

    # annotate gene names, remove noncoding if required. 
    df['gene'] = df['info'].apply(get_gene)
    if not allow_noncoding:
        noncoding = set(df[df['type']!='gene']['gene'].unique())
        df = df[~df['gene'].isin(noncoding)].copy()
    df = df[df['type'].isin(['gene', 'pseudogene'])].copy()
    
    df['span'] = df['end'] - df['start']
    df = df[['chr', 'start', 'end', 'span', 'gene']].copy()
    return df 

def load_gff_gencode(filepath: str, allow_noncoding: bool, allow_x: bool=False, allow_y: bool=False) -> pd.DataFrame:

    def get_name(text: str) -> str:
        PATTERN = r'^.*?gene_name=([A-Za-z0-9_\.-]+).*?$'
        m = re.match(PATTERN, text)
        assert m is not None
        return m.group(1)

    def get_gtype(text: str) -> str:
        PATTERN = r'^.*?gene_type=([A-Za-z0-9_-]+).*?$'
        m = re.match(PATTERN, text)
        assert m is not None
        return m.group(1)

    data = []
    with open(filepath, 'r') as fp:
        line = fp.readline().strip('\n')
        i = 0
        while line:
            if line.startswith('#'):
                line = fp.readline().strip('\n')
                continue 
            
            i += 1
            if i % 10_000 == 0:
                print(f"processed {i} records...", end='\r')
            
            lsplit = line.split('\t')
            chrom = lsplit[0].replace('chr', '')

            if lsplit[2] != 'gene':
                line = fp.readline().strip('\n')
                continue
            if chrom == 'X' and not allow_x:
                line = fp.readline().strip('\n')
                continue 
            if chrom == 'Y' and not allow_y:
                line = fp.readline().strip('\n')
                continue 
            
            start = int(lsplit[3])
            end = int(lsplit[4])
            gene = get_name(lsplit[8])
            gtype = get_gtype(lsplit[8])

            if gtype == 'protein_coding' or allow_noncoding:
                data.append((chrom, start, end, gene))

            line = fp.readline().strip('\n')

    print(f"processed {i} records... done")
    df = pd.DataFrame.from_records(data, columns=['chr', 'start', 'end', 'gene'])
    df['chr'] = df['chr'].astype(str)
    df['start'] = df['start'].astype(int)
    df['end'] = df['end'].astype(int)
    df['gene'] = df['gene'].astype(str)
    return df



############################
### CNA EVENT VALIDATORS ###
############################

class CNAEventValidator(ABC):

    @abstractmethod
    def validate(self, seg: pd.Series, wgd: bool) -> bool:
        ...
    
    def is_loh(self, seg: pd.Series) -> bool:
        if seg.chr == 'X':
            return False
        return seg.nMaj >= 1 and seg.nMin == 0

class AmpValidator(CNAEventValidator):
    
    def validate(self, seg: pd.Series, wgd: bool) -> bool:
        loh = self.is_loh(seg)

        if seg.chr == 'X':
            return seg.tcn >= 3
        elif loh and wgd:
            return seg.tcn >= 6
        elif wgd:
            return seg.tcn >= 7
        else:
            return seg.tcn >= 5

class ShallowDelValidator(CNAEventValidator):
    
    def validate(self, seg: pd.Series, wgd: bool) -> bool:

        if seg.chr == 'X':
            return seg.tcn == 0
        elif wgd:
            return seg.tcn <= 1
        else:
            return seg.tcn == 0

class DeepDelValidator(CNAEventValidator):
    
    def validate(self, seg: pd.Series, wgd: bool) -> bool:
        return seg.tcn == 0
    

##################
### EXTRACTION ###
##################

def extract_refseq(segs: pd.DataFrame, gff: pd.DataFrame, args: argparse.Namespace) -> pd.DataFrame:
    validator_lut = dict()
    if args.allow_amp: 
        validator_lut['CNA↑'] = AmpValidator
    if args.allow_shallow_del: 
        validator_lut['CNA↓'] = ShallowDelValidator
    if args.allow_deep_del:
        validator_lut['CNA↓↓'] = DeepDelValidator

    data = []
    # iterate allowed variant types
    for vtype, validator_c in validator_lut.items():
        # get segments which meet CNA variant criteria
        df = segs.copy()
        validator = validator_c()
        df['PASS'] = df.apply(validator.validate, wgd=args.wgd, axis=1)
        df = df[df['PASS']==True].copy()

        # get genes which intersect valid segments & update data
        for idx, row in df.iterrows():
            coords = row['chr'] + ':' + str(row['start']) + '-' + str(row['end'])
            genes = get_genes(row, gff, args.min_span)
            for gene in genes:
                data.append((coords, '.', '.', gene, vtype, 'CNAgene', round(row['frac_ccf'], 2)))

    table = pd.DataFrame.from_records(data, columns=['coords', 'ref', 'alt', 'gene', 'vtype', 'annotation', 'est_ccf'])
    return table

def get_genes(seg: pd.Series, gff: pd.DataFrame, min_span: float) -> set[str]:
    gslice = gff[(gff['chr']==seg.chr) & (gff['start']<seg.end) & (gff['end']>seg.start)]
    out = set()
    for rec in gslice.itertuples():
        thresh = min_span * rec.span
        overlap = min(rec.end, seg.end) - max(rec.start, seg.start)
        if overlap >= thresh:
            out.add(rec.gene) 
    return out

def extract_gencode(segs: pd.DataFrame, gff: pd.DataFrame, args: argparse.Namespace) -> pd.DataFrame:
    shared_chroms = sorted(list(set(segs['chr'].unique()) & set(gff['chr'].unique())))
    
    validator_lut = dict()
    if args.allow_amp: 
        validator_lut['CNA↑'] = AmpValidator
    if args.allow_shallow_del: 
        validator_lut['CNA↓'] = ShallowDelValidator
    if args.allow_deep_del:
        validator_lut['CNA↓↓'] = DeepDelValidator

    data = []
    # iterate allowed variant types
    for vtype, validator_c in validator_lut.items():
        validator = validator_c()

        # iter chroms (performance)        
        for chrom in shared_chroms:
            gslice = gff[gff['chr']==chrom]
            sslice = segs[segs['chr']==chrom]

            # get segments in this chromosome which are valid (meets the CNA criteria)
            sslice['PASS'] = sslice.apply(validator.validate, wgd=args.wgd, axis=1)
            sslice = sslice[sslice['PASS']==True]
            
            for gene in gslice.itertuples():
                cnas = sslice[(sslice['start']<gene.end) & (sslice['end']>gene.start)].copy()

                if cnas.shape[0] == 0:
                    continue 
                
                if chrom == 'X':
                    print()
                    print(cnas)
                cnas['start'] = cnas['start'].clip(lower=gene.start)
                cnas['end'] = cnas['end'].clip(upper=gene.end)
                cnas['span'] = cnas['end'] - cnas['start']
                gspan = gene.end - gene.start
                sspan = cnas['span'].sum()
                thresh = args.min_span * gspan
                if sspan >= thresh:
                    coords = gene.chr + ':' + str(gene.start) + '-' + str(gene.end)
                    ccf = round(cnas['frac_ccf'].mean(), 2)
                    data.append((coords, '.', '.', gene.gene, vtype, 'CNAgene', ccf))

    table = pd.DataFrame.from_records(data, columns=['coords', 'ref', 'alt', 'gene', 'vtype', 'annotation', 'est_ccf'])
    return table


def load_cmdline_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description='Extracts CNA events from a PPCG donor SCNA file.')
    parser.add_argument('--scna', type=str, required=True, help='Path to Battenberg SCNA segments (for CCF estimation).')
    parser.add_argument('--gff', type=str, required=True, help='Path to GFF file containing gene coordinates.')
    parser.add_argument('--loh', action='store_true', help='Whether to extract LOH events.')
    parser.add_argument('--wgd', action='store_true', help='Whether the sample has clonal WGD.')
    parser.add_argument('--allow-subclonal', action='store_true', help='Include segments where two states are observed (subclonal).')
    parser.add_argument('--allow-amp', action='store_true', help='Include amplifications.')
    parser.add_argument('--allow-shallow-del', action='store_true', help='Include shallow deletions.')
    parser.add_argument('--allow-deep-del', action='store_true', help='Include deep deletions.')
    parser.add_argument('--allow-x', action='store_true', help='Include chromosome X variants.')
    parser.add_argument('--allow-noncoding', action='store_true', help='Include non-coding variants.')
    parser.add_argument('--min-span', type=float, default=0.99, help='The minimum proportion of the gene/peak affected to call a CNA variant.')
    parser.add_argument('--outfile', type=str, required=True, help='Path to output file containing extracted variants')
    args = parser.parse_args()
    return args

if __name__ == '__main__':
    main()


# def load_cmdline_args()
    # parser.add_argument('--cosmic', type=str, required=True, help='Path to COSMIC cancer genes file.')
    # parser.add_argument('--peaks', type=str, required=True, help='Path to recurrent Gistic2 peaks file.')
    # parser.add_argument('--max-qval', type=float, default=0.05, help='The maximum Gistic2 Q-value for a peak to be included.')

# def main()
    # old approach using select regions (cosmic genes / gistic peaks)
    # gene_ivals = gene_ivals.drop_duplicates('label').copy()
    # peak_ivals = load_CNApeak_intervals(args)
    # gene_ivals = load_CNAgene_intervals(args)
    # gframe['vtype'] = 'CNAgene'
    # pframe['vtype'] = 'CNApeak'
    # df = pd.concat([gframe, pframe], ignore_index=True)

# def load_CNAgene_intervals(args: argparse.Namespace) -> pd.DataFrame:
#     cframe = _load_cosmic_intervals(args)
#     pframe = _load_peak_intervals(args)

#     df = pd.DataFrame()
#     for rec in pframe.itertuples():
#         cslice = cframe[cframe['chr']==rec.chr].copy()
#         cslice = cslice[(cslice['start']<rec.end) & (cslice['end']>rec.start)]
#         cslice['allowed'] = rec.allowed
#         df = pd.concat([df, cslice], ignore_index=True)
#     return df 

# def load_CNApeak_intervals(args: argparse.Namespace) -> pd.DataFrame:
#     return _load_peak_intervals(args)
    
# def _load_cosmic_intervals(args: argparse.Namespace) -> pd.DataFrame:
#     df = pd.read_csv(args.cosmic, sep='\t', header=0)
#     df = df[df['GENOME_START'].notna()]
#     df['CHROMOSOME'] = df['CHROMOSOME'].astype(str).apply(lambda x: f"chr{x}")
#     df['GENOME_START'] = df['GENOME_START'].astype(int)
#     df['GENOME_STOP'] = df['GENOME_STOP'].astype(int)
#     df = df[['CHROMOSOME', 'GENOME_START', 'GENOME_STOP', 'GENE_SYMBOL']]
#     df = df.rename(columns={
#         'CHROMOSOME': 'chr',
#         'GENOME_START': 'start',
#         'GENOME_STOP': 'end',
#         'GENE_SYMBOL': 'label',
#     })
#     df['chr'] = df['chr'].apply(lambda x: x.replace('chr', ''))
#     return df.copy()

# def _load_peak_intervals(args: argparse.Namespace) -> pd.DataFrame:
#     df = pd.read_csv(args.peaks, sep='\t', header=0)
#     df['chr'] = df['Wide_Peak_Limits'].apply(lambda x: x.split(':')[0].replace('chr', ''))
#     df['start'] = df['Wide_Peak_Limits'].apply(lambda x: x.split(':')[1].split('-')[0])
#     df['end'] = df['Wide_Peak_Limits'].apply(lambda x: x.split(':')[1].split('-')[1])
#     df['start'] = df['start'].astype(int)
#     df['end'] = df['end'].astype(int)
#     df = df[df['qvalues']<=args.max_qval].copy()
#     df = df.rename(columns={'Variant_Classification': 'allowed', 'Cytoband': 'label'})
#     df = df[['chr', 'start', 'end', 'label', 'allowed']]
#     return df.copy()



# def extract_event(segs: pd.DataFrame, etype: str, q_chrom: str, q_start: int, q_end: int, wgd: bool, min_span: float):
#     """
#     detects whether the donor has this type of CN event within the query interval. 
#     all available samples are checked. 
#     """
#     df = segs.copy()

#     # intersecting segments
#     df = df[df['chr']==q_chrom]
#     df = df[(df['start']<q_end) & (df['end']>q_start)].copy()

#     # clipping segments to interval
#     df['start'] = df['start'].apply(lambda x: max(x, q_start))
#     df['end'] = df['end'].apply(lambda x: min(x, q_end))

#     # filter segments to only those which qualify for this event
#     vsegs = validate_segments(df, etype, wgd)
#     vsegs = vsegs[vsegs['PASS']==True].copy()

#     # early exit if cumulative length of filtered segments too short
#     vsegs['span'] = vsegs['end'] - vsegs['start']
#     cumspan = vsegs['span'].sum()
#     thresh = int((q_end-q_start) * min_span)
#     if cumspan <= thresh:
#         return False, None, 0.0
    
#     coords = f"{q_chrom}:{q_start}-{q_end}"
#     maxspan = vsegs['span'].max()
#     ccf = vsegs[vsegs['span']==maxspan]['frac_ccf'].values[0]
#     return True, coords, ccf

# def validate_segments(segs: pd.DataFrame, etype: str, wgd: bool) -> pd.DataFrame:
#     amap = {
#         'AMP': AmpValidator,
#         'DEL': ShallowDelValidator,
#         'DEEP DEL': DeepDelValidator,
#     }
#     validator_c = amap[etype]
#     validator = validator_c()
#     segs['PASS'] = segs.apply(validator.validate, wgd=wgd, axis=1)
#     return segs