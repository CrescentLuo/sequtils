import argparse
from utils import reverse_complement, check_intersect, one_based_start
import twobitreader
import swifter
import os
import pandas as pd
from tqdm import tqdm
import cgranges as cr


def get_args():
    parser = argparse.ArgumentParser(
        description='Compute the impact of novel exons')
    parser.add_argument(
        '-i', '--intron', required=True,
        help='input intron annotation')
    parser.add_argument(
        '-b', '--bed', required=True,
        help='table of novel exons in BED format')
    parser.add_argument(
        '-o', '--output', default='./insertion_intron_idx.tab', help='output file')
    args = parser.parse_args()
    return args


def build_cgranges(df):
    g = cr.cgranges()
    max_intron_dict = dict()
    for idx, row in tqdm(df.iterrows(), total=df.shape[0]):
        seqname = row['seqname']
        start = row['start']
        end = row['end']
        name = row['name']
        ts_id = name.split('_')[0]
        intron_num = int(name.split('_')[-1])
        g.add(seqname, start, end, idx)
        if ts_id not in max_intron_dict:
            max_intron_dict[ts_id] = 0
        if intron_num > max_intron_dict[ts_id]:
            max_intron_dict[ts_id] = intron_num
    g.index()
    return g,max_intron_dict

def calculate_intron_len(intron):
    intron['length'] = intron['end'] - intron['start']
    return intron

if __name__ == '__main__':
    args = get_args()
    gtf_intron = pd.read_csv(
        args.intron,
        sep='\t',names=['seqname', 'start', 'end', 'name', 'score', 'strand'],
        dtype={'start': int, 'end': int})
    gtf_intron = gtf_intron.swifter.apply(calculate_intron_len)
    # cr_ts = build_cgranges(gtf_ts)
    cr_intron,max_intron_dict = build_cgranges(gtf_intron)
    
    # convert to dict
    gtf_intron = gtf_intron.to_dict('index')

    exon_target = pd.read_csv(
        args.bed, sep='\t',
        names=['chrom', 'start', 'end', 'name', 'score', 'strand'])
    exon_target['intron_idx'] =-1
    exon_target['intron_idx_rev'] =-1
    exon_target['first/last intron'] = 'False'    
    exon_target['intron_len'] = -1            
    for idx, exon in tqdm(exon_target.iterrows(), total=exon_target.shape[0]):
        # for idx, exon in exon_target.iterrows():
        chrom = exon['chrom']
        start = exon['start']
        end = exon['end']
        pos = 0
        for intron_s, intron_end, intron_idx in cr_intron.overlap(chrom, start, end):
            pos = intron_idx
        intron_name = gtf_intron[pos]['name']
        ts_id = intron_name.split('_')[0]
        intron_idx = int(intron_name.split('_')[1])
        intron_idx_rev = max_intron_dict[ts_id] - intron_idx + 1
        exon_target.at[idx,'intron_idx'] = intron_idx
        exon_target.at[idx,'intron_idx_rev'] = intron_idx_rev
        exon_target.at[idx,'intron_len'] = gtf_intron[pos]['intron_len']
        if intron_idx == 1:
            exon_target.at[idx,'first/last intron'] = 'first'
        if intron_idx == 2:
            exon_target.at[idx,'first/last intron'] = 'second'
        if intron_idx_rev == 1:
            exon_target.at[idx,'first/last intron'] = 'last'
        if intron_idx_rev == 2:
            exon_target.at[idx,'first/last intron'] = 'second to last'
        

    print(exon_target.head())
    exon_target.to_csv(args.output, sep='\t', index=False)
