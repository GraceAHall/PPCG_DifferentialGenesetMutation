
import pandas as pd 


def load_scna(filepath: str, allow_subclonal: bool) -> pd.DataFrame:
    df = pd.read_csv(filepath, sep='\t', header=0)
    df = df.rename(columns={'startpos': 'start', 'endpos': 'end'})
    df['chr'] = df['chr'].astype(str)
    df['start'] = df['start'].astype(int)
    df['end'] = df['end'].astype(int)
    df = df[['chr', 'start', 'end', 'nMaj1_A', 'nMin1_A', 'frac1_A', 'nMaj2_A', 'nMin2_A', 'frac2_A']]
    
    data = []
    for idx, rec in df.iterrows():
        segments = []
        if not rec.isna()['nMaj2_A']:
            majority = max(rec.frac1_A, rec.frac2_A)
        else:
            majority = rec.frac1_A
        segments.append((
            rec.chr, 
            rec.start, 
            rec.end, 
            rec.nMaj1_A, 
            rec.nMin1_A, 
            rec.nMaj1_A+rec.nMin1_A,
            rec.frac1_A / majority
        ))
        if not rec.isna()['nMaj2_A']:
            segments.append((
                rec.chr, 
                rec.start, 
                rec.end, 
                rec.nMaj2_A, 
                rec.nMin2_A, 
                rec.nMaj2_A+rec.nMin2_A,
                rec.frac2_A / majority
            ))
        if len(segments)==2 and not allow_subclonal:
            continue 
        else:
            data += segments

    cnframe = pd.DataFrame.from_records(data, columns=['chr', 'start', 'end', 'nMaj', 'nMin', 'tcn', 'frac_ccf'])
    return cnframe


class CCFestimator:
    def __init__(self, purity: float, scna_path: str) -> None:
        self.purity = purity
        self.cnframe = load_scna(scna_path, allow_subclonal=True)

    def est_ccf(self, chrom: str, position: int, vaf: float) -> float:
        f = vaf 
        p = self.purity

        # get copy number information
        cnslice = self.cnframe[self.cnframe['chr']==chrom].copy()
        cnslice = cnslice[(cnslice['start']>=position) & (cnslice['end']<=position)].copy()
        if cnslice.shape[0] == 0:
            # assume normal copy number. 
            # tcn == 2 for autosomes, 1 for sex chromosomes
            segs = [[2.0, 1.0]] if chrom.isdigit() else [[1.0, 1.0]]
        else:
            segs = [[float(rec.tcn), rec.frac_ccf] for rec in cnslice.itertuples()] # type: ignore

        data = []
        for nt, frac_ccf in segs:
            m = max(1, round(( f/p ) * ( p*nt + 2*(1-p) )))
            ccf = ( f/(m*p) ) * ( p*nt + 2*(1-p) )
            diff = abs(frac_ccf-ccf)
            data.append((round(ccf, 2), diff))
        data.sort(key=lambda x: x[-1])
        return min(1, data[0][-2])
    
