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
    for idx, row in tqdm(df.iterrows(), total=df.shape[0]):
        seqname = row['seqname']
        start = row['start']
        end = row['end']
        g.add(seqname, start, end, idx)
    g.index()
    return g


if __name__ == '__main__':
    args = get_args()
    gtf_intron = pd.read_csv(
        args.intron,
        sep='\t',names=['seqname', 'start', 'end', 'name', 'score', 'strand'],
        dtype={'start': int, 'end': int})
    # cr_ts = build_cgranges(gtf_ts)
    cr_intron = build_cgranges(gtf_intron)
    
    # convert to dict
    gtf_intron = gtf_intron.to_dict('index')

    exon_target = pd.read_csv(
        args.bed, sep='\t',
        names=['chrom', 'start', 'end', 'name', 'score', 'strand'])
    exon['intron_idx'] =-1
    for idx, exon in tqdm(exon_target.iterrows(), total=exon_target.shape[0]):
        # for idx, exon in exon_target.iterrows():
        chrom = exon['chrom']
        start = exon['start']
        end = exon['end']
        for intron_s, intron_end, intron_idx in cr_intron.overlap(chrom, start, end):
            break
        intron_name = gtf_intron[intron_idx]['name']
        intron_idx = int(intron_name.split('_')[-1])
        exon['intron_idx'] = intron_idx
    exon_target.to_csv(args.output, sep='\t', index=False)
