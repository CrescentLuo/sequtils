import pandas as pd
import swifter
import itertools
import re

# check chrom

def check_chrom_re():
    ref_chrom_re = '^chr([1-9]|[12][0-2]|[X]|[Y])$'
    return ref_chrom_re

def one_based_start(s):
    return s+1

def zero_based_start(s):
    return s-1
# nucleotides sequence utils

def reverse_complement(seq):
    # return reverse complementary sequence of input
    try:
        assert set(seq).issubset(set("ACGTN"))
    except AssertionError:
        raise ValueError(
            "Sequence {} contains invalid values: {}"
            .format(seq, set(seq) - set('ACGTN'))
        )
    return ''.join('TGCAN'['ACGTN'.index(s)] for s in seq.upper()[::-1])

# Annotation Parser
#   Function: parse_gtf
#   Input - GTF file from GENCODE annotation database
#   Output:
#       1. Gene-level features: gene_id, gene_type, gene_name extracted from attributes column
#       2. Transcript-level features: gene_id, gene_type, gene_name, transcript_id, 
#           transcript_support_level extracted from attributes column, cds_start and cds_end 
#           assigned the same as transcript_start.
#       3. Exon-level features: gene_id, gene_type, gene_name, transcript_id, exon_id,
#           exon_number extracted from attributes column.
#       4. CDS-level features: same as exon-level features.

def parse_gtf(gtf_file):
    gtf = pd.read_csv(
        gtf_file, sep='\t', comment='#',
        names=['seqname', 'source', 'feature', 'start', 'end', 'score',
               'strand', 'frame', 'attribute'])
    gtf_gene = gtf[gtf['feature'] == 'gene'].swifter.apply(
        parse_gtf_attrs, axis=1)
    gtf_transcript = gtf[gtf['feature'] ==
                         'transcript'].swifter.apply(parse_gtf_attrs, axis=1)
    gtf_exon = gtf[gtf['feature'] == 'exon'].swifter.apply(
        parse_gtf_attrs, axis=1)
    gtf_CDS = gtf[gtf['feature'] == 'CDS'].swifter.apply(
        parse_gtf_attrs, axis=1)
    return gtf_gene, gtf_transcript, gtf_exon, gtf_CDS


def parse_gtf_attrs(feature):
    feature_type = feature['feature']
    attrs = feature['attribute']
    feature_start = feature['start']
    regex = r"(\S+\s\S+;)"
    attrs_dict = dict()
    for attr in re.finditer(regex, attrs):
        attr_key = attr[1].split()[0]
        attr_val = attr[1].split()[1].rstrip(';').strip('"')
        attrs_dict[attr_key] = attr_val

    # parse attributes from gene features
    # gene_id, gene_type, gene_name
    if feature_type == 'gene':
        feature['gene_id'] = attrs_dict['gene_id']
        feature['gene_type'] = attrs_dict['gene_type']
        feature['gene_name'] = attrs_dict['gene_name']
    # parse attributes from transcript features
    # gene_id, gene_type, gene_name
    if feature_type == 'transcript':
        feature['gene_id'] = attrs_dict['gene_id']
        feature['gene_type'] = attrs_dict['gene_type']
        feature['gene_name'] = attrs_dict['gene_name']
        feature['cds_start'] = feature_start
        feature['cds_end'] = feature_start
        feature['transcript_id'] = attrs_dict['transcript_id']
        if 'transcript_support_level' in attrs_dict:
            feature['transcript_support_level'] = attrs_dict['transcript_support_level']
        else:
            feature['transcript_support_level'] = 'NA'
    if feature_type == 'exon':
        feature['gene_id'] = attrs_dict['gene_id']
        feature['gene_type'] = attrs_dict['gene_type']
        feature['gene_name'] = attrs_dict['gene_name']
        feature['transcript_id'] = attrs_dict['transcript_id']
        feature['exon_id'] = attrs_dict['exon_id']
        feature['exon_number'] = int(attrs_dict['exon_number'])
    if feature_type == 'CDS':
        feature['gene_id'] = attrs_dict['gene_id']
        feature['gene_type'] = attrs_dict['gene_type']
        feature['gene_name'] = attrs_dict['gene_name']
        feature['transcript_id'] = attrs_dict['transcript_id']
        feature['exon_id'] = attrs_dict['exon_id']
        feature['exon_number'] = int(attrs_dict['exon_number'])

    return feature

# todo:Extract introns from transcript 
# Extract alternative exons from annotation
# input gtf_exon is an instance of pandas dataframe
# input gtf_ts 
def get_alt_exons(gtf_ts, gtf_exon):
    gtf_exon = gtf_exon.to_dict('index')
    gtf_ts = gtf_ts.to_dict('index')
    gene_ts_cnt = dict()
    exon_usage = dict()
    ts_exon_dict = dict()
    for idx in gtf_exon:
        transcript_id = gtf_exon[idx]['transcript_id']
        if transcript_id not in ts_exon_dict:
            ts_exon_dict[transcript_id] = list()
        ts_exon_dict[transcript_id].append(idx)
    
    for idx in gtf_ts:
        gene_id = gtf_ts[idx]['gene_id']
        if gene_id not in gene_ts_cnt:
            gene_ts_cnt[gene_id] = 0
        gene_ts_cnt[gene_id] += 1
    
    for transcript_id in ts_exon_dict:
        exon_set = sorted(ts_exon_dict[transcript_id])
        exon_num = len(exon_set)
        for i, exon_idx in enumerate(exon_set):
            if i == 0 or i == exon_num - 1:
                continue

            gene_id = gtf_exon[exon_idx]['gene_id']
            exon_id = gtf_exon[exon_idx]['seqname'] + str(gtf_exon[exon_idx]['start']) + str(gtf_exon[exon_idx]['end']) + gtf_exon[exon_idx]['strand']
            if exon_id not in exon_usage:
                exon_usage[exon_id] = 0
            exon_usage[exon_id] += 1
    alt_exon_idx = list()
    alt_exon_coord = set()
    for idx in gtf_exon:
        exon_id = gtf_exon[idx]['seqname'] + str(gtf_exon[idx]['start']) + str(gtf_exon[idx]['end']) + gtf_exon[idx]['strand']
        gene_id = gtf_exon[idx]['gene_id']
        if exon_id in exon_usage and gene_id in gene_ts_cnt:
            if exon_usage[exon_id] != gene_ts_cnt[gene_id]:
                if exon_id not in alt_exon_coord:
                    alt_exon_idx.append(idx)
                    alt_exon_coord.add(exon_id)
                
    return alt_exon_idx


# extract introns from annotation
def get_introns(gtf_ts, gtf_exon):
    gtf_ts = gtf_ts.to_dict('index')
    gtf_exon = gtf_exon.to_dict('index')
    ts_exon_dict = dict()
    for idx in gtf_exon:
        transcript_id = gtf_exon[idx]['transcript_id']
        if transcript_id not in ts_exon_dict:
            ts_exon_dict[transcript_id] = list()
        ts_exon_dict[transcript_id].append(idx)
    
    for ts_idx in gtf_ts:
        transcript_id = gtf_ts[ts_idx]['transcript_id']
        exon_set = sorted(ts_exon_dict[transcript_id])
        exon_num = len(exon_set)
        if exon_num < 2:
            continue
        intron_chrom = gtf_ts[ts_idx]['seqname']
        intron_strand = gtf_ts[ts_idx]['strand']
        intron_start = 0
        intron_end = 0
        for i, exon_idx in enumerate(exon_set):
            if i == 0:
                if intron_strand == '+':
                    intron_start = gtf_exon[exon_idx]['end']
                else:
                    intron_end = gtf_exon[exon_idx]['start']
            else:
                if intron_strand == '+':
                    intron_end = gtf_exon[exon_idx]['start']
                else:
                    intron_start = gtf_exon[exon_idx]['end']
                intron_name = transcript_id+'_'+str(i)
                yield '{}\t{}\t{}\t{}\t{}\t{}\n'.format(intron_chrom, intron_start, intron_end, intron_name, '.', intron_strand)
                if intron_strand == '+':
                    intron_start = gtf_exon[exon_idx]['end']
                else:
                    intron_end = gtf_exon[exon_idx]['start']

# Calculate the correct CDS region at trascript level.
def assign_cds_to_ts(gtf_ts, gtf_cds):
    gtf_ts = gtf_ts.to_dict('index')
    gtf_cds = gtf_cds.to_dict('index')
    ts_cds_dict = dict()

    for idx in gtf_cds:
        transcript_id = gtf_cds[idx]['transcript_id']
        if transcript_id not in ts_cds_dict:
            ts_cds_dict[transcript_id] = list()
        ts_cds_dict[transcript_id].append(idx)
    for idx in gtf_ts:
        transcript_id = gtf_ts[idx]['transcript_id']
        if transcript_id in ts_cds_dict:
            for cds_idx in ts_cds_dict[transcript_id]:
                cds_start = int(gtf_cds[cds_idx]['start'])
                cds_end = int(gtf_cds[cds_idx]['end'])
                ts_cds_start = int(gtf_ts[idx]['cds_start'])
                ts_cds_end = int(gtf_ts[idx]['cds_end'])
                ts_start = int(gtf_ts[idx]['start'])
                if ts_cds_start == ts_start:
                    gtf_ts[idx]['cds_start'] = cds_start
                    gtf_ts[idx]['cds_end'] = cds_end
                else:
                    if cds_start < ts_cds_start:
                        gtf_ts[idx]['cds_start'] = cds_start
                    if cds_end > ts_cds_end:
                        gtf_ts[idx]['cds_end'] = cds_end
    gtf_ts = pd.DataFrame.from_dict(gtf_ts, orient='index')
    return gtf_ts

# Coordinates transformation
def check_intersect(coord_1, coord_2, zero_based_1=True, zero_based_2=True):
    if zero_based_1:
        coord_1[0] += 1
    if zero_based_2:
        coord_2[0] += 1
    if coord_1[1] < coord_2[0] or coord_1[0] > coord_2[1]:
        return False
    else:
        return True

def build_ts_dict(df):
    ts_dict = dict()
    for idx in df:
        transcript_id = df[idx]['transcript_id']
        if transcript_id not in ts_dict:
            ts_dict[transcript_id] = list()
        ts_dict[transcript_id].append(idx)
    return ts_dict

def concat_interval_seq(ts_id, ts_dict, interval_set, genome_seq):
    concat_seq = ''
    strand = '+'
    for interval_id in ts_dict[ts_id]:
        chrom = interval_set[interval_id]['seqname']
        if interval_set[interval_id]['strand'] == '-':
            strand = '-'
        start = int(interval_set[interval_id]['start'])
        start = zero_based_start(start)
        end = int(interval_set[interval_id]['end'])
        concat_seq = concat_seq + genome_seq.fetch(chrom,start,end)
    if strand == '-':
        return reverse_complement(concat_seq)
    return concat_seq
    
def ambiguity_nt(seq):
    '''
    N = A or C or G or T (any)
    B = C or G or T (not A)
    D = A or G or T (not C)
    H = A or C or T (not G)
    V = A or C or G (not T)
    W = A or T (weak)
    S = C or G (strong)
    R = A or G (purine)
    Y = C or T (pyrimidine)
    M = A or C (amino)
    K = G or T (keto)
    '''
    ambg_dict = {
        'N':['A', 'C', 'G', 'T'],
        'B':['C', 'G', 'T'],
        'D':['A', 'G', 'T'],
        'H':['A', 'C', 'T'],
        'V':['A', 'C', 'G'],
        'W':['A', 'T'],
        'S':['C', 'G'],
        'R':['A', 'G'],
        'Y':['C', 'T'],
        'M':['A', 'C'],
        'K':['G', 'T'],
        'A':['A'],
        'T':['T'],
        'C':['C'],
        'G':['G']
    }
    combinations = list()
    for nt in seq:
        combinations.append(ambg_dict[nt])
    seq_permutations = list()
    for id, comb_seq in enumerate(itertools.product(*combinations)):
        comb_seq = ''.join(comb_seq)
        seq_permutations.append(comb_seq)
    return seq_permutations