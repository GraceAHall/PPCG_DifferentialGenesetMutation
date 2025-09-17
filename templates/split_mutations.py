
import argparse 
import pandas as pd 


def main() -> None:
    args = load_cmdline_args()

    muts = pd.read_csv(args.mutations, sep='\t', header=0)
    sheet = pd.read_csv(args.samplesheet, sep='\t', header=0)
    pdonors = set(sheet[sheet['cohort']==args.posclass]['donor'].unique())
    posmuts = muts[muts['donor'].isin(pdonors)].copy()
    negmuts = muts[~muts['donor'].isin(pdonors)].copy()
    posmuts.to_csv(args.outfile_pos, sep='\t', index=False)
    negmuts.to_csv(args.outfile_neg, sep='\t', index=False)


def load_cmdline_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description='Splits mutations into two seperate files based on contrast.')
    parser.add_argument('--posclass', type=str, required=True, help='Which group in samplesheet is considered the positive class.')
    parser.add_argument('--mutations', type=str, required=True, help='Path to merged mutations.')
    parser.add_argument('--samplesheet', type=str, required=True, help='Path to samplesheet mapping donors to groups.')
    parser.add_argument('--outfile-pos', type=str, required=True, help='Output file for positive class mutations.')
    parser.add_argument('--outfile-neg', type=str, required=True, help='Output file for negative class mutations.')
    args = parser.parse_args()
    return args

if __name__ == '__main__':
    main()