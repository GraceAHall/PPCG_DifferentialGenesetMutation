
import argparse 
import pandas as pd 
from utils import filter_hypermutators

def main() -> None:
    args = load_cmdline_args()
    df = pd.read_csv(args.mutations, sep='\t', header=0)
    df = df.drop('sample', axis=1)
    
    # hypermutators
    df_filt, results = filter_hypermutators(df, args.zscore_thresh)
    results = results.reset_index()

    # write to file.
    df_filt.to_csv(args.outfile_mutations, sep='\t', index=False)
    results.to_csv(args.outfile_summary, sep='\t', index=False)

def load_cmdline_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description='Merges all variant types into single file for a given PPCG donor.')
    parser.add_argument('--mutations', type=str, required=True, help='Path to input file containing merged variants (tsv).')
    parser.add_argument('--zscore-thresh', type=float, required=True, help='outlier threshold via Modified Z-Score (MAD) (standard deviations from the mean).')
    parser.add_argument('--outfile-mutations', type=str, required=True, help='Path to output .tsv file containing filtered variants.')
    parser.add_argument('--outfile-summary', type=str, required=True, help='Path to output .tsv file containing hypermutator summary.')
    args = parser.parse_args()
    return args

if __name__ == '__main__':
    main()
