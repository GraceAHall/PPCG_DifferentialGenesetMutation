
import re 
import numpy as np
import pandas as pd 
import argparse


def main() -> None:
    args = load_cmdline_args()
    df = load_svs(args.vcf)
    sframe = merge_pairs(df)
    sframe.to_csv(args.outfile, sep='\t', index=False)

def merge_pairs(df: pd.DataFrame) -> pd.DataFrame:
    data = []
    df = df.sort_values('ID')
    for pair in df['PAIR'].unique():
        dfslice = df[df['PAIR']==pair]
        chr1 = dfslice[dfslice['ID']==f"{pair}_1"]['CHROM'].values[0]
        pos1 = dfslice[dfslice['ID']==f"{pair}_1"]['POS'].values[0]
        chr2 = dfslice[dfslice['ID']==f"{pair}_2"]['CHROM'].values[0]
        pos2 = dfslice[dfslice['ID']==f"{pair}_2"]['POS'].values[0]
        coords = f"{chr1}:{pos1}-{chr2}:{pos2}"
        vaf = dfslice['VAF'].values[0]
        dfslice = dfslice.drop_duplicates(subset=['OUTCOME', 'GENES'])
        for rec in dfslice.itertuples():
            data.append((rec.ID, coords, rec.GENES, rec.SVCLASS, rec.OUTCOME, vaf, np.nan))
    
    FIELDS = ['id', 'coords', 'gene', 'vtype', 'annotation', 'VAF', 'est_ccf']
    vframe = pd.DataFrame.from_records(data=data, columns=FIELDS)
    return vframe

def load_svs(filepath: str) -> pd.DataFrame:
    snpeff_format = ['Allele', 'Annotation', 'Annotation_Impact', 'Gene_Name', 'Gene_ID', 'Feature_Type', 'Feature_ID', 'Transcript_BioType', 'Rank', 'HGVS.c', 'HGVS.p', 'cDNA.pos / cDNA.length', 'CDS.pos / CDS.length', 'AA.pos / AA.length', 'Distance', 'ERRORS / WARNINGS / INFO']
    SV_SNPEFF_FORMAT_MAP = {f: i for i, f in enumerate(snpeff_format)}
    
    def _get_snpeff_ann_values(info: str, field: str) -> list[str]:
        pattern = r'.*ANN=(.*)'
        m = re.match(pattern, info)
        if m is None:
            return []
        annotations = [x.split('|') for x in m.group(1).split(',')]
        idx = SV_SNPEFF_FORMAT_MAP[field]
        values = [x[idx] for x in annotations]
        values = ['.' if x=='' else x for x in values] 
        return values

    def _get_brass_delly_svclass(info: str) -> str:
        pattern = r'.*;SVCLASS=(\w+);.*'
        m = re.match(pattern, info)
        assert m is not None 
        svclass = m.group(1)
        if svclass in ['h2hINV', 't2tINV']:
            svclass = 'INV'
        return svclass
    
    # how many lines to skip (headers)
    with open(filepath, 'r') as fp:
        lines = fp.readlines()
        for i, ln in enumerate(lines):
            if not ln.startswith('##'):
                n_skip = i
                break

    # load vcf
    table = pd.read_csv(filepath, sep='\t', skiprows=n_skip+1, header=None)
    table.columns = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', 'NORMAL', 'TUMOUR']
    table = table[table['FILTER']=='PASS']

    data = []
    for idx, row in table.iterrows():
        chrom = row['CHROM']
        pos = row['POS']
        ident = row['ID']
        info = row['INFO']

        svclass = _get_brass_delly_svclass(info)
        vaf = np.nan

        all_outcomes = _get_snpeff_ann_values(info, 'Annotation')
        all_impacts = _get_snpeff_ann_values(info, 'Annotation_Impact')
        all_genes = _get_snpeff_ann_values(info, 'Gene_Name')
        assert len(all_outcomes) == len(all_impacts)
        assert len(all_impacts) == len(all_genes)
        for outcome, impact, genes in zip(all_outcomes, all_impacts, all_genes):
            # removing IGKC, weird bug????
            if 'IGKC' in genes:
                continue 
            if impact != 'HIGH':
                continue 
            data.append((chrom, pos, ident, svclass, vaf, outcome, genes))
    df = pd.DataFrame.from_records(data=data, columns=['CHROM', 'POS', 'ID', 'SVCLASS', 'VAF', 'OUTCOME', 'GENES'])
    df['PAIR'] = df['ID'].apply(lambda x: x.split('_')[0])
    return df

def load_cmdline_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description='Extracts SVs from a PPCG donor VCF file.')
    parser.add_argument('--vcf', type=str, required=True, help='Path to structural variant VCF file.')
    parser.add_argument('--outfile', type=str, required=True, help='Path to output file containing extracted variants')
    args = parser.parse_args()
    return args

if __name__ == '__main__':
    main()

