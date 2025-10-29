
import argparse 
import pandas as pd 

def main() -> None:
    args = load_cmdline_args()
    df = pd.read_csv(args.mutations, sep='\t', header=0)
    df = df.drop('sample', axis=1)
    df.to_csv(args.outfile, sep='\t', index=False)

def load_cmdline_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description='Merges all variant types into single file for a given PPCG donor.')
    parser.add_argument('--mutations', type=str, required=True, help='Path to input file containing merged variants (tsv).')
    parser.add_argument('--outfile', type=str, required=True, help='Path to output file containing filtered variants.')
    args = parser.parse_args()
    return args

if __name__ == '__main__':
    main()
