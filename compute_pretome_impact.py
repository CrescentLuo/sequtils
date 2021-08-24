import argparse
import bisect
import twobitreader
from operator import itemgetter, attrgetter
import swifter
import os
import pandas as pd
from tqdm import tqdm
import cgranges as cr
import dask


def get_args():
    parser = argparse.ArgumentParser(
        description='Compute the impact of novel exons')
    parser.add_argument(
        '-d', '--gtfdir', required=True,
        help='input annotation files in GTF format')
    parser.add_argument(
        '-b', '--bed', required=True,
        help='table of novel exons in BED format')
    parser.add_argument('--genome', required=False,
                        help='genome sequence file in 2bit format')
    args = parser.parse_args()
    return args


def search_stop_codon(seq, lo=0):
    # stop codon TAG TAA TGA
    stop_codon = set("TAG", "TAA", "TGA")
    for i in range(0, len(seq), 3):
        codon = seq[i:i+2]
        if codon in stop_codon:
            return i, codon
    return -1, "NA"


def track_ts(exon, cr_ts, gtf_ts, gtf_exon, gtf_ts_exon, gtf_ts_cds, gtf_cds):
    exon_s = exon['start']
    exon_e = exon['end']
    exon_chrom = exon['chrom']
    exon_strand = exon['strand']
    # get the corresponding gene
    for ts_s, ts_e, ts_idx in cr_ts.overlap(
            exon_chrom, exon_s, exon_e):
        transcript_id = gtf_ts[ts_idx]['transcript_id']
        ts_strand = gtf_ts[ts_idx]['strand']
        if exon_strand != ts_strand:
            continue
        # 检测target 是否在coding region
        
        c1_idx = -1
        c2_idx = -1
        last_eej_idx = -1
        for idx in gtf_ts_exon[transcript_id]:
            ts_exon_s = gtf_exon[idx]['start']
            ts_exon_e = gtf_exon[idx]['end']
            #print("{}\t{}<--{}-->{}\t{}".format(ts_exon_s,ts_exon_e,exon_strand,exon_s,exon_e))
            if exon_strand == '+':
                last_eej_idx = idx
            else:
                if last_eej_idx == -1:
                    last_eej_idx =idx
            if ts_exon_e < exon_s:
                if exon_strand == '+':
                    c1_idx = idx
                else:
                    c2_idx = idx
            if ts_exon_s > exon_e:
                if exon_strand == '+':
                    c2_idx = idx
                else:
                    c1_idx = idx
        if c1_idx != -1 and c2_idx !=-1:

            # 
    return -1, -1


def get_dis_last_eej():
    return


def get_cds():
    return


def search_ptc():
    return


def build_cgranges(df):
    g = cr.cgranges()
    for idx, row in tqdm(df.iterrows(), total=df.shape[0]):
        seqname = row['seqname']
        start = row['start']
        end = row['end']
        g.add(seqname, start, end, idx)
    g.index()
    return g

def build_ts_exon_dict(df):
    ts_exon_dict = dict()
    
    for idx in df:
        transcript_id = df[idx]['transcript_id']
        if transcript_id not in ts_exon_dict:
            ts_exon_dict[transcript_id] = list()
        ts_exon_dict[transcript_id].append(idx)
    return ts_exon_dict

if __name__ == '__main__':
    args = get_args()
    print(args.gtfdir)
    gtf_ts = pd.read_csv(
        os.path.join(
            args.gtfdir, 'gencode.v31.primary_assembly.annotation.ts.tab'),
        sep='\t', dtype={'start': int, 'end': int})
    gtf_ts = gtf_ts[gtf_ts['gene_type'] == 'protein_coding']
    gtf_ts = gtf_ts[gtf_ts['transcript_support_level'].isin(['1',
                                                             '2'])]
    cr_ts = build_cgranges(gtf_ts)
    gtf_ts = gtf_ts.to_dict('index')

    gtf_exon = pd.read_csv(
        os.path.join(
            args.gtfdir, 'gencode.v31.primary_assembly.annotation.exon.tab'),
        sep='\t', dtype={'start': int, 'end': int})
    gtf_exon = gtf_exon.to_dict('index')
    gtf_ts_exon = build_ts_exon_dict(gtf_exon)
    gtf_cds = pd.read_csv(
        os.path.join(
            args.gtfdir, 'gencode.v31.primary_assembly.annotation.CDS.tab'),
        sep='\t', dtype={'start': int, 'end': int})
    gtf_cds = gtf_cds.to_dict('index')
    gtf_ts_cds = build_ts_exon_dict(gtf_cds)
    
    #genome = twobitreader.TwoBitFile(args.genome)
    exon_target = pd.read_csv(
        args.bed, sep='\t',
        names=['chrom', 'start', 'end', 'name', 'score', 'strand'])
    for idx, exon in tqdm(exon_target.iterrows(), total=exon_target.shape[0]):
    #for idx, exon in exon_target.iterrows():
        c1_idx, c2_idx = track_ts(exon, cr_ts=cr_ts, gtf_ts=gtf_ts, gtf_exon=gtf_exon, gtf_ts_exon=gtf_ts_exon, gtf_ts_cd=gtf_ts_cds, gtf_cds=gtf_cds)
        if c1_idx != -1 and c2_idx != -1:
            print("##################")
            print(gtf_exon[c1_idx])
            print(exon)
            print(gtf_exon[c2_idx])
