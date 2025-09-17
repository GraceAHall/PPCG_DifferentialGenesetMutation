
import re 
import argparse
import pandas as pd 
from collections import defaultdict
# from utils import expand_fusions


def main() -> None:
    args = load_cmdline_args()
    muts = pd.read_csv(args.mutations, sep='\t', header=0)
    sizes = pd.read_csv(args.sizes, sep='\t', header=0)
    gsets = pd.read_csv(args.genesets, sep='\t', header=0)
    gsets = gsets.rename(columns={'Gene': 'gene', 'Pathway': 'geneset'})

    muts = muts[muts['vtype']!='CNApeak'].copy()

    sizes = standardise(sizes, muts, args)
    gsets = standardise(gsets, muts, args)
    sizes = standardise(sizes, gsets, args)
    gsets = gsets[['gene', 'geneset']].copy()
    sizes = sizes.drop_duplicates('gene')
    # print(sizes[sizes.duplicated('gene', keep=False)])

    # ensure all mutated genes have entry in sizes
    # ensure all geneset genes have entry in sizes
    mut_genes = sorted(list(muts['gene'].unique()))
    gset_genes = sorted(list(gsets['gene'].unique()))
    sizes = update_sizes(sizes, mut_genes)
    sizes = update_sizes(sizes, gset_genes)

    sizes.to_csv(args.outfile_sizes, sep='\t', index=False)
    gsets.to_csv(args.outfile_gsets, sep='\t', index=False)


PATTERN_INT = r'(\w+)-IT\d+'
PATTERN_DVG = r'(\w+)-DT'
PATTERN_ASS = r'([\w-]+)-AS\d'
PATTERN_NUM = r'(\w+)-\d+'
PATTERN_THRU = r'(\w+)-(\w\w\w+)'

MIR_SIZE  = 22
SNOR_SIZE = 150
RNU_SIZE  = 150

def update_sizes(df: pd.DataFrame, genelist: list[str]) -> pd.DataFrame:
    med_exon = int(df['cum_exon_len'].median())
    med_cds = int(df['cum_cds_len'].median())
    med_span = int(df['span'].median())
    # print(df.head())
    # print()
    # print(df.tail())
    # print()
    # print(df[df.duplicated('gene', keep=False)])

    df = df.set_index('gene')
    data = []
    for gene in genelist:
        size_genes = set(df.index.to_list())
        if gene in size_genes:
            continue
        if gene.startswith('MT-'):
            continue 

        m_dvg = re.match(PATTERN_DVG, gene)
        m_int = re.match(PATTERN_INT, gene)
        m_ass = re.match(PATTERN_ASS, gene)
        m_num = re.match(PATTERN_NUM, gene)
        
        # divergent transcripts
        if m_dvg is not None:
            base = m_dvg.group(1)
            if base in size_genes:
                data.append((
                    df.loc[base, 'cum_exon_len'], 
                    df.loc[base, 'cum_cds_len'], 
                    df.loc[base, 'span'], 
                    gene, 
                ))
                continue 
        
        # intronic transcripts
        if m_int is not None:
            base = m_int.group(1)
            if base in size_genes:
                data.append((
                    df.loc[base, 'cum_exon_len'], 
                    df.loc[base, 'cum_cds_len'], 
                    df.loc[base, 'span'], 
                    gene, 
                ))
                continue 
        
        # antisense transcripts
        if m_ass is not None:
            base = m_ass.group(1)
            if base in size_genes:
                data.append((
                    df.loc[base, 'cum_exon_len'], 
                    df.loc[base, 'cum_cds_len'], 
                    df.loc[base, 'span'], 
                    gene, 
                ))
                continue 
        
        # gene subunit numbers
        if m_num is not None:
            base = m_num.group(1)
            if base in size_genes:
                data.append((
                    df.loc[base, 'cum_exon_len'], 
                    df.loc[base, 'cum_cds_len'], 
                    df.loc[base, 'span'], 
                    gene, 
                ))
                continue 
        
        # gene subunit numbers
        if m_num is not None:
            base = m_num.group(1)
            if base in size_genes:
                data.append((
                    df.loc[base, 'cum_exon_len'], 
                    df.loc[base, 'cum_cds_len'], 
                    df.loc[base, 'span'], 
                    gene, 
                ))
                continue
        
        # miRNA: 22
        if gene.startswith('MIR'):
            data.append((MIR_SIZE, MIR_SIZE, MIR_SIZE, gene))
            continue 
        if gene.startswith('SNOR'):
            data.append((SNOR_SIZE, SNOR_SIZE, SNOR_SIZE, gene))
            continue 
        if gene.startswith('RNU'):
            data.append((RNU_SIZE, RNU_SIZE, RNU_SIZE, gene))
            continue 
        
        # fallback: average lengths
        data.append((med_exon, med_cds, med_span, gene))

    df = df.reset_index('gene')
    df = df[['cum_exon_len', 'cum_cds_len', 'span', 'gene']].copy()
    df_new = pd.DataFrame.from_records(data, columns=['cum_exon_len', 'cum_cds_len', 'span', 'gene'])
    sizes = pd.concat([df, df_new], ignore_index=True)
    return sizes


def split_genebases(genes: set[str]) -> set[str]:
    out = set()
    for gene in genes:
        bases = get_genebases(gene)
        out.update(bases)
    return out

def get_genebases(gene: str) -> list[str]:
    m_dvg = re.match(PATTERN_DVG, gene)
    m_int = re.match(PATTERN_INT, gene)
    m_num = re.match(PATTERN_NUM, gene)
    m_ass = re.match(PATTERN_ASS, gene)
    m_thru = re.match(PATTERN_THRU, gene)

    if m_num is not None:
        return [gene]
    elif m_int is not None:
        return [gene]
    elif m_dvg is not None:
        # print(f"diverging: {gene}")
        return [m_dvg.group(1)]
    elif m_ass is not None:
        # print(f"antisense: {gene}")
        return [m_ass.group(1)]
    elif m_thru is not None:
        # print(f"readthrough: {gene}")
        return gene.split('-')
    else:
        return [gene]

def join_genebases(src: set[str], dest: set[str], lut: dict[str, str]) -> dict[str, str]:
    out = dict()
    for src_gene in src:
        # shared
        if src_gene in dest:
            out[src_gene] = src_gene
            continue
            
        # mapped
        if src_gene in lut:
            out[src_gene] = lut[src_gene]
            continue 

        # base mapped?
        bases = get_genebases(src_gene)
        if all([x in lut for x in bases]):
            new_sym_base = '-'.join([lut[x] for x in bases])
            
            if len(bases) > 1:
                # try readthrough
                if new_sym_base in dest:
                    print(f'readthrough: {src_gene} -> {new_sym_base}')
                    out[src_gene] = new_sym_base
                    continue

            # try divering
            new_sym_dt = f"{new_sym_base}-DT"
            if new_sym_dt in dest:
                print(f'diverging: {src_gene} -> {new_sym_dt}')
                out[src_gene] = new_sym_dt
                continue

            # try antisense
            for i in range(10):
                new_sym_ass = f"{new_sym_base}-AS{i}"
                if new_sym_ass in dest:
                    print(f'antisense: {src_gene} -> {new_sym_ass}')
                    out[src_gene] = new_sym_ass
                    break 

    return out

def standardise(df_src: pd.DataFrame, df_dest: pd.DataFrame, args: argparse.Namespace) -> pd.DataFrame:
    ss = SymbolStandardiser(args.gff, args.gtf, args.hgnc)
    ss.load()

    src_genes = set(df_src['gene'].unique())
    dest_genes = set(df_dest['gene'].unique())

    src_bases = split_genebases(src_genes)
    dest_bases = split_genebases(dest_genes)

    lut = ss.map2(src_bases, dest_bases)
    mapping = join_genebases(src_genes, dest_genes, lut)

    df = df_src.copy()
    df['mapped'] = df['gene'].apply(lambda x: mapping.get(x, x))
    assert df['mapped'].isna().sum() == 0
    df = df.drop('gene', axis=1)
    df = df.rename(columns={'mapped': 'gene'})
    return df
        


class SymbolStandardiser:

    _map_cache: dict[str, set[str]]
    
    _official_symbols: set[str]
    _official_ids: set[str]

    _hgnc_refid2sym: dict[str, str]
    _hgnc_id2sym: dict[str, str]
    _hgnc_sym2id: dict[str, str]

    _hgnc_prev2curr: dict[str, set[str]]
    _hgnc_alt2curr: dict[str, set[str]]

    _bio_refsym_2_refid: dict[str, set[str]]
    _gff_refsym_2_refid: dict[str, set[str]]
    
    _bio_refsym_2_hgncid: dict[str, set[str]]
    _gff_refsym_2_hgncid: dict[str, set[str]]

    def __init__(self, gff_path: str, gtf_path: str, hgnc_path: str) -> None:
        self.gff_path = gff_path
        self.gtf_path = gtf_path
        self.hgnc_path = hgnc_path

    def load(self) -> None:
        print('Loading symbol standardisation LUTS... ', end='\r')
        # self._map_cache: dict[str, set[str]] = dict()
        self._load_refseq_gff()
        self._load_biomart_gtf()
        self._load_hgnc()
        print('Loading symbol standardisation LUTS... done')

    ### PUBLIC METHODS ###
    def map2(self, src: set[str], dest: set[str]) -> dict[str, str]:
        mapping: dict[str, str] = dict()
        shared = src & dest
        src_uniq = src - dest
        dest_uniq = dest - src

        # standardise dest genes once to avoid quadratic complexity
        backmap: dict[str, set[str]] = defaultdict(set)
        for dest_sym in dest_uniq:
            dest_current = list(self.current_symbols(dest_sym))
            for dest_curr in dest_current:
                if dest_curr != dest_sym:
                    backmap[dest_curr].add(dest_sym)
    
        # no mapping needed
        for sym in shared:
            mapping[sym] = sym

        # requies mapping
        for src_sym in src_uniq:
            res = self._do_map(src_sym, backmap, dest_uniq)
            if res is not None:
                mapping[src_sym] = res

        # print()
        # print(f"src symbols:  {len(src)}") 
        # print(f"src unique:  {len(src_uniq)}") 
        # print(f"dest symbols:  {len(dest)}") 
        # print(f"dest unique: {len(dest_uniq)}") 
        # print(f"shared:  {len(shared)}") 
        # print(f"total mapped: {len([k for k, v in mapping.items() if v is not None])}") 
        # print(f"total unmapped: {len([k for k, v in mapping.items() if v is None])}") 
        # print(f"src coverage:   {len(src_covered)/len(src)*100:.1f}%")
        # print(f"dest coverage:   {len(dest_covered)/len(dest)*100:.1f}%")
        # print(f"identical symbols: {identical}") 
        # print(f"different symbols: {different}") 

        return mapping

    def _do_map(self, query: str, backmap: dict, ref_uniq: set) -> str|None:
            
        # dest is outdated / alias
        if query in backmap:
            return ','.join(sorted(list(backmap[query])))
        
        # map src to current
        query_current = self.current_symbols(query)

        # src has single current symbol 
        if len(query_current) == 1:
            curr = query_current.pop()
            
            # src is outdated / alias 
            if curr in ref_uniq:
                return curr

            # src and dest are outdated / alias 
            if curr in backmap:
                return ','.join(sorted(list(backmap[curr])))

        # src has multiple current symbols
        if len(query_current) > 1:
            shared = query_current & ref_uniq
            if len(shared) == 1:
                return shared.pop()
            elif len(shared) > 1:
                return ','.join(sorted(list(shared)))
                # print(query, query_current, shared)
                # raise NotImplementedError

            bmhits = {}
            for curr in query_current:
                if curr in backmap:
                    bmhits[curr] = backmap[curr]
            if len(bmhits) == 1:
                dest_symbols = bmhits.popitem()[1]
                return ','.join(sorted(list(dest_symbols)))
            elif len(bmhits) > 1:
                print(query, query_current, bmhits)
                raise NotImplementedError

        # assume src is simply not in dest symbols.
        return None
    
    def current_symbols(self, query: str) -> set[str]:
        # if query in self._map_cache:
        #     return self._map_cache[query]
        
        current = set()

        # active HGNC
        if query in self._official_symbols:
            current = set([query])

        # using HGNC previous symbol or alias symbol
        if not current:
            current = self._get_hgnc(query)
        
        # using BioMART symbol -> refSeq ID -> HGNC ID
        if not current:
            current = self._get_biomart(query)
        
        # using refSeq symbol -> refSeq ID -> HGNC ID
        if not current:
            current = self._get_refseq(query)

        # self._map_cache[query] = current
        return current
    

    ### PRIVATE METHODS ###
    def _get_hgnc(self, query: str) -> set[str]:
        if query in self._hgnc_prev2curr:
            return self._hgnc_prev2curr[query]
        elif query in self._hgnc_alt2curr:
            return self._hgnc_alt2curr[query]
        return set()

    def _get_biomart(self, query: str) -> set[str]:
        hgnc_ids = self._bio_refsym_2_hgncid[query]
        ref_ids = self._bio_refsym_2_refid[query]

        out = set()
        if hgnc_ids:
            for uid in hgnc_ids:
                if uid in self._hgnc_id2sym:
                    out.add(self._hgnc_id2sym[uid])
        elif ref_ids:
            for uid in ref_ids:
                if uid in self._hgnc_refid2sym:
                    out.add(self._hgnc_refid2sym[uid])
        
        return set()
    
    def _get_refseq(self, query: str) -> set[str]:
        hgnc_ids = self._gff_refsym_2_hgncid[query]
        ref_ids = self._gff_refsym_2_refid[query]

        out = set()
        if hgnc_ids:
            for uid in hgnc_ids:
                if uid in self._hgnc_id2sym:
                    out.add(self._hgnc_id2sym[uid])
        elif ref_ids:
            for uid in ref_ids:
                if uid in self._hgnc_refid2sym:
                    out.add(self._hgnc_refid2sym[uid])
        
        return set()

    def _load_biomart_gtf(self) -> None:
        self._bio_refsym_2_refid = defaultdict(set) 
        self._bio_refsym_2_hgncid = defaultdict(set) 
        
        mart = pd.read_csv(self.gtf_path, sep='\t', header=0, dtype='str')
        for i, row in mart.iterrows():
            ref_sym = row['Gene name']
            ref_id = row['NCBI gene (formerly Entrezgene) ID']
            hgnc_id = row['HGNC ID']
            if isinstance(ref_id, str):
                self._bio_refsym_2_refid[ref_sym].add(ref_id)
            if isinstance(hgnc_id, str):
                self._bio_refsym_2_hgncid[ref_sym].add(hgnc_id)

    def _load_refseq_gff(self) -> None:
        self._gff_refsym_2_refid = defaultdict(set)
        self._gff_refsym_2_hgncid = defaultdict(set)
        with open(self.gff_path, 'r') as fp:
            line = fp.readline().strip()
            while line:
                if not line.startswith('#'):
                    ref_id = None
                    ref_sym = None
                    hgnc_id = None
                    lsplit = line.split('\t') 
                    ftype = lsplit[2]
                    if ftype == 'gene':
                        details = lsplit[8]
                        for item in details.split(';'):
                            if item.startswith('gene='):
                                ref_sym = item.split('=')[-1]
                            elif item.startswith('Dbxref='):
                                for subitem in item.replace('Dbxref=', '').split(','):
                                    if subitem.startswith('HGNC'):
                                        hgnc_id = subitem.split(':')[-1]
                                    elif subitem.startswith('GeneID'):
                                        ref_id = subitem.split(':')[-1]
                        
                        assert ref_sym is not None
                        assert ref_id is not None
                        self._gff_refsym_2_refid[ref_sym].add(ref_id)
                        if hgnc_id is not None:
                            self._gff_refsym_2_hgncid[ref_sym].add(hgnc_id)
                line = fp.readline().strip()

    def _load_hgnc(self) -> None:
        hgnc = pd.read_csv(self.hgnc_path, sep='\t', quotechar='"', header=0, dtype=str)
        hgnc['hgnc_id'] = hgnc['hgnc_id'].apply(lambda x: x.split(':', 1)[1])

        # list of active hgnc symbols
        self._official_symbols = set(hgnc['symbol'].unique())

        # list of active hgnc ids
        self._official_ids = set(hgnc['hgnc_id'].unique())

        # ids <-> symbols
        self._hgnc_refid2sym = dict()
        self._hgnc_id2sym = dict()
        self._hgnc_sym2id = dict()
        for i, row in hgnc.iterrows():
            sym = row['symbol']
            self._hgnc_refid2sym[row['entrez_id']] = sym
            self._hgnc_id2sym[row['hgnc_id']] = sym
            self._hgnc_sym2id[sym] = row['hgnc_id']

        # lookup tables for gene symbol aliases/changes.
        self._hgnc_prev2curr = self._load_prevsymbol_lut(hgnc)
        self._hgnc_alt2curr = self._load_altsymbol_lut(hgnc)
    
    def _load_prevsymbol_lut(self, hgnc: pd.DataFrame) -> dict[str, set[str]]:
        lut = defaultdict(set)
        
        # previous names
        for i, row in hgnc[hgnc['prev_symbol'].notna()].iterrows():
            current = row['symbol']
            alts = row['prev_symbol'].split('|')
            assert len(alts) >= 1
            for alt in alts:
                lut[alt].add(current)
        
        return lut
    
    def _load_altsymbol_lut(self, hgnc: pd.DataFrame) -> dict[str, set[str]]:
        lut = defaultdict(set)
        
        # aliases
        for i, row in hgnc[hgnc['alias_symbol'].notna()].iterrows():
            current = row['symbol']
            alts = row['alias_symbol'].split('|')
            assert len(alts) >= 1
            for alt in alts:
                lut[alt].add(current)
        
        return lut


def load_cmdline_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description='Standardises symbols.')
    parser.add_argument('--mutations', type=str, required=True, help='Path to extracted mutations file (tsv).')
    parser.add_argument('--genesets', type=str, required=True, help='Path to GSEA genesets file (tsv).')
    parser.add_argument('--sizes', type=str, required=True, help='Path to gene sizes (tsv).')
    parser.add_argument('--gff', type=str, required=True, help='Path to refseq GRCh37 file (gff).')
    parser.add_argument('--gtf', type=str, required=True, help='Path to ensembl GRCh37 file (gtf or txt).')
    parser.add_argument('--hgnc', type=str, required=True, help='Path to HGNC file (txt).')
    parser.add_argument('--outfile-sizes', type=str, required=True, help='Path to output file mapping genes to sizes.')
    parser.add_argument('--outfile-gsets', type=str, required=True, help='Path to output file mapping genes to genesets.')

    args = parser.parse_args()
    return args

if __name__ == '__main__':
    main()