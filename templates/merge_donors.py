
import re 
import sys 
import argparse
import pandas as pd 
from glob import glob 

pd.set_option('display.max_rows', 1000)
pd.set_option('display.max_columns', 1000)
pd.set_option('display.width', 1000)
pd.set_option('display.max_colwidth', 1000)
# pd.options.display.float_format = '{:.2f}'.format

def main() -> None:
    args = load_cmdline_args()
    
    # init the merged table
    merged = merge_files(args)
    merged_exp = expand_fusions(merged)

    # report summary stats
    log_summary(merged, merged_exp, args)

    # save to file
    merged_exp.to_csv(args.outfile, index=False, sep='\t')

def expand_fusions(df: pd.DataFrame) -> pd.DataFrame:
    exp = df.copy()
    exp['gene'] = exp['gene'].apply(lambda x: x.split('&'))
    exp = exp.explode('gene')
    return exp.copy()

def log_summary(df: pd.DataFrame, df_exp: pd.DataFrame, args: argparse.Namespace) -> None:
    original_stdout = sys.stdout

    # open logfile, write all print statements to file (stdout capture)
    basename = args.outfile.rsplit('.', 1)[0]
    with open(f"{basename}.log", 'w') as f:
        sys.stdout = f

        print('------------------')
        print('--- Basic Info ---')
        print('------------------')
        print(f"{df['donor'].nunique()} donors")
        print(f"{df['gene'].nunique()} genes")
        print(f"{df['coords'].nunique()} variants")
        
        print()
        print('-----------------------------------')
        print('--- Gene Summary (top 50 genes) ---')
        print('-----------------------------------')
        
        ndonors = df['donor'].nunique()
        counts = pd.DataFrame(index=sorted(list(df['gene'].unique())))
        counts['variants'] = df['gene'].value_counts()
        counts['donors'] = df.groupby('gene')['donor'].nunique()
        counts['donors prop.'] = counts['donors'] / ndonors
        counts['donors prop.'] = counts['donors prop.'].apply(lambda x: f"{x*100:.1f} %")
        counts = counts.sort_values('donors', ascending=False)
        print()
        print(counts.head(50))

        print()
        print('-----------------------------------------------------')
        print('--- Gene Summary (top 50 genes, fusions expanded) ---')
        print('-----------------------------------------------------')

        ndonors = df_exp['donor'].nunique()
        counts = pd.DataFrame(index=sorted(list(df_exp['gene'].unique())))
        counts['variants'] = df_exp['gene'].value_counts()
        counts['donors'] = df_exp.groupby('gene')['donor'].nunique()
        counts['donors prop.'] = counts['donors'] / ndonors
        counts['donors prop.'] = counts['donors prop.'].apply(lambda x: f"{x*100:.1f} %")
        counts = counts.sort_values('donors', ascending=False)
        print()
        print(counts.head(50))

        print()
        print('-----------------------------')
        print('--- Variant Class Summary ---')
        print('-----------------------------')
        
        ndonors = df['donor'].nunique()
        counts = pd.DataFrame(index=sorted(list(df['vclass'].unique())))
        counts['variants'] = df['vclass'].value_counts()
        counts['donors'] = df.groupby('vclass')['donor'].nunique()
        counts['donors prop.'] = counts['donors'] / ndonors
        counts['donors prop.'] = counts['donors prop.'].apply(lambda x: f"{x*100:.1f} %")
        counts = counts.sort_values('donors', ascending=False)
        print()
        print(counts)
        
        print()
        print('-----------------------------')
        print('--- Variant Type Summary ---')
        print('-----------------------------')
        
        ndonors = df['donor'].nunique()
        labels = ['Sequence Variants', 'Structural Variants', 'Copy Number']
        vclasses = [['SNV', 'INDEL'], ['SV'], ['CNA']]
        for label, vclass_l in zip(labels, vclasses):
            dfslice = df[df['vclass'].isin(vclass_l)].copy()
            counts = pd.DataFrame(index=sorted(list(dfslice['vtype'].unique())))
            counts['variants'] = dfslice['vtype'].value_counts()
            counts['donors'] = dfslice.groupby('vtype')['donor'].nunique()
            counts['donors prop.'] = counts['donors'] / ndonors
            counts['donors prop.'] = counts['donors prop.'].apply(lambda x: f"{x*100:.1f} %")
            counts = counts.sort_values('donors', ascending=False)
            print()
            print(f"({label})")
            print(counts)

        print()
        print('----------------------------------')
        print('--- Variant Annotation Summary ---')
        print('----------------------------------')
 
        ndonors = df['donor'].nunique()
        for vclass in sorted(list(df['vclass'].unique())):
            dfslice = df[df['vclass']==vclass].copy()
            counts = pd.DataFrame(index=sorted(list(dfslice['annotation'].unique())))
            counts['variants'] = dfslice['annotation'].value_counts()
            counts['donors'] = dfslice.groupby('annotation')['donor'].nunique()
            counts['donors prop.'] = counts['donors'] / ndonors
            counts['donors prop.'] = counts['donors prop.'].apply(lambda x: f"{x*100:.1f} %")
            counts = counts.sort_values('donors', ascending=False)
            print()
            print(f"({vclass})")
            print(counts)

        print()
        print('### NOTE: FROM HERE FUSIONS HAVE BEEN EXPANDED ###')
        print('### Eg "ERG&TMPRSS2" will be expanded to individual records for "ERG" and "TMPRSS2" ###')
        print()

        print('--------------------------------------------------')
        print('--- Genes (top 30) Stratified by Variant Class ---')
        print('--------------------------------------------------')
        print('Numbers represent donors.')

        temp = df_exp.drop_duplicates(['gene', 'donor', 'vclass']).copy()
        counts = temp.groupby('gene')['vclass'].value_counts().unstack().fillna(0).astype(int)
        counts['donors'] = temp.groupby('gene')['donor'].nunique()
        counts = counts.sort_values('donors', ascending=False)
        print()
        print(counts.head(30))
        
        print()
        print('-------------------------------------------------------------------')
        print('--- Genes (top 20 per variant class) Stratified by Variant Type ---')
        print('-------------------------------------------------------------------')
        print('Numbers represent donors.')

        labels = ['Sequence Variants', 'Structural Variants', 'Copy Number']
        vclasses = [['SNV', 'INDEL'], ['SV'], ['CNA']]
        for label, vclass_l in zip(labels, vclasses):
            dfslice = df_exp[df_exp['vclass'].isin(vclass_l)].copy()
            dfslice = dfslice.drop_duplicates(['gene', 'donor']).copy()
            genes = dfslice['gene'].value_counts().head(20).index.to_list()
            dfslice = dfslice[dfslice['gene'].isin(genes)]
            counts = dfslice.groupby('gene')['vtype'].value_counts().unstack().fillna(0).astype(int)
            counts['donors'] = dfslice.groupby('gene')['donor'].nunique()
            counts = counts.sort_values('donors', ascending=False)
            print()
            print(f"({label})")
            print(counts)

        print()
        print('-------------------------------------------------------------------------')
        print('--- Genes (top 20 per variant class) Stratified by Variant Annotation ---')
        print('-------------------------------------------------------------------------')
        print('Numbers represent donors.')

        for vclass in sorted(list(df['vclass'].unique())):
            dfslice = df_exp[df_exp['vclass']==vclass].copy()
            dfslice = dfslice.drop_duplicates(['gene', 'donor']).copy()
            genes = dfslice['gene'].value_counts().head(20).index.to_list()
            dfslice = dfslice[dfslice['gene'].isin(genes)]
            counts = dfslice.groupby('gene')['annotation'].value_counts().unstack().fillna(0).astype(int)
            counts['donors'] = dfslice.groupby('gene')['donor'].nunique()
            counts = counts.sort_values('donors', ascending=False)
            print()
            print(f"({vclass})")
            print(counts)

    sys.stdout = original_stdout

def merge_files(args: argparse.Namespace) -> pd.DataFrame:
    FIELDS = ['donor', 'coords', 'ref', 'alt', 'gene', 'vclass', 'vtype', 'annotation', 'est_ccf']
    merged = pd.DataFrame(columns=FIELDS)

    # include SVs for all donors
    for filepath in glob(f"{args.svdir}/*.tsv"):
        donor = filepath.split('/')[-1].split('.')[0][:8]
        df = pd.read_csv(filepath, header=0, sep='\t')
        df['donor'] = donor
        df['vclass'] = 'SV'
        df['ref'] = '.'
        df['alt'] = '.'
        df = df[FIELDS].copy()
        merged = pd.concat([merged, df], ignore_index=True)
    
    # include CNA events for all donors
    for filepath in glob(f"{args.cnadir}/*.tsv"):
        donor = filepath.split('/')[-1].split('.')[0][:8]
        df = pd.read_csv(filepath, header=0, sep='\t')
        df['donor'] = donor
        df['vclass'] = 'CNA'
        # df = df.rename(columns={'msubtype': 'annotation', 'mtype': 'vtype'})
        df = df[FIELDS].copy()
        merged = pd.concat([merged, df], ignore_index=True)
    
    # include SNVs for all donors
    for filepath in glob(f"{args.snvdir}/*.tsv"):
        donor = filepath.split('/')[-1].split('.')[0][:8]
        df = pd.read_csv(filepath, header=0, sep='\t')
        df['donor'] = donor
        df['vclass'] = 'SNV'
        df['vtype'] = 'SNV'
        df = df.rename(columns={'REF': 'ref', 'ALT': 'alt'})
        df = df[FIELDS].copy()
        merged = pd.concat([merged, df], ignore_index=True)
    
    # include INDELs for all donors
    for filepath in glob(f"{args.indeldir}/*.tsv"):
        donor = filepath.split('/')[-1].split('.')[0][:8]
        df = pd.read_csv(filepath, header=0, sep='\t')
        df['donor'] = donor
        df['vclass'] = 'INDEL'
        df['vtype'] = 'INDEL'
        df = df.rename(columns={'REF': 'ref', 'ALT': 'alt'})
        df = df[FIELDS].copy()
        merged = pd.concat([merged, df], ignore_index=True)

    return merged


def load_cmdline_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description='Merges all variant types into single file for a given PPCG donor.')
    parser.add_argument('--svdir', type=str, required=True, help='Path to directory containing extracted structural variants (tsv).')
    parser.add_argument('--snvdir', type=str, required=True, help='Path to directory containing extracted SNVs (tsv).')
    parser.add_argument('--cnadir', type=str, required=True, help='Path to directory containing extracted CNA events (tsv).')
    parser.add_argument('--indeldir', type=str, required=True, help='Path to directory containing extracted INDELs (tsv).')
    parser.add_argument('--outfile', type=str, required=True, help='Path to output file containing merged variants.')
    parser.add_argument('--logfile', type=str, required=True, help='Path to output log file containing summary information.')
    args = parser.parse_args()
    return args

if __name__ == '__main__':
    main()

