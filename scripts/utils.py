import pandas as pd
import textwrap
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

def write_seq(seq_set, file_path):
    with open(file_path, 'w') as outfile:
        for seq_id in seq_set:
            outfile.write(seq_id+'\n')
            outfile.write(textwrap.fill(seq_set[seq_id], width=80) + '\n')