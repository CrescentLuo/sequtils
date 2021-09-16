import argparse
import bisect
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
        '-d', '--gtfdir', required=True,
        help='input annotation files in GTF format')
    parser.add_argument(
        '-b', '--bed', required=True,
        help='table of novel exons in BED format')
    parser.add_argument('--genome', required=False,
                        help='genome sequence file in 2bit format')
    parser.add_argument(
        '-o', '--output', default='./proteome_impact.tab', help='output file')
    args = parser.parse_args()
    return args


def search_ptc(seq, lo=0):
    # stop codon TAG TAA TGA
    stop_codon = set(["TAG", "TAA", "TGA"])
    for i in range(0, len(seq), 3):
        codon = seq[i:i+3]
        if codon in stop_codon:
            return i, codon
    return -1, -1

# extract genomic sequence based on coordinate, can be both 0-based or 1-based
def extract_seq(genome, chrom, start, end, zero_based=True):
    if not zero_based:
        start -= 1
    return genome[chrom][start:end]


def get_residue(cds, cds_frame):
    coden_len = 3
    strand = cds['strand']
    return (cds['end'] - cds['start'] - cds_frame + 1) % coden_len

def check_frame_shift(exon_seq):
    coden_len = 3
    if len(exon_seq) % coden_len != 0:
        #return('frame shift_{}_{}'.format(len(exon_seq, transcript_id)))
        return('ORF-disrupting')
    else:
        return('ORF-preserving)')


def track_ts(
    exon, genome, cr_ts, gtf_ts, gtf_exon, gtf_ts_exon, gtf_ts_cds,
        gtf_cds):
    coden_len = 3
    exon_s = int(exon['start'])
    exon_e = int(exon['end'])
    exon_chrom = exon['chrom']
    exon_strand = exon['strand']
    # get the corresponding gene
    impact_isoform = set()
    
    for ts_s, ts_e, ts_idx in cr_ts.overlap(
            exon_chrom, exon_s, exon_e):
        transcript_id = gtf_ts[ts_idx]['transcript_id']
        ts_strand = gtf_ts[ts_idx]['strand']
        # require the exon and ts from the same strand
        if exon_strand != ts_strand:
            continue
        # skip transcripts without cds region (cds_start == cds_end)
        if gtf_ts[ts_idx]['cds_start'] == gtf_ts[ts_idx]['cds_end']:
            impact_isoform.add('ORF-preserving)')
            continue
        elif not check_intersect([exon_s, exon_e], [gtf_ts[ts_idx]['cds_start'],gtf_ts[ts_idx]['cds_end']], 1,0):
            impact_isoform.add('ORF-preserving)')
            continue
        c1_idx = -1
        c2_idx = -1
        last_eej_idx = -1
        for idx in sorted(gtf_ts_cds[transcript_id]):
            ts_exon_s = gtf_cds[idx]['start']
            ts_exon_e = gtf_cds[idx]['end']
            if exon_strand == '+':
                last_eej_idx = idx
            else:
                if last_eej_idx == -1:
                    last_eej_idx = idx
            if ts_exon_e < one_based_start(exon_s):
                if exon_strand == '+':
                    c1_idx = idx
                elif c1_idx==-1:
                    c1_idx = idx
            if ts_exon_s > exon_e:
                if exon_strand == '-':
                    c2_idx = idx
                elif c2_idx==-1:
                    c2_idx = idx
                    
        if c1_idx != -1 and c2_idx != -1:
            c1_cds = gtf_cds[c1_idx]
            c2_cds = gtf_cds[c2_idx]
            exon_seq = ""
            c2_frame = int(c2_cds['frame'])
            c1_frame = int(c1_cds['frame'])
            c1_residue = get_residue(c1_cds, c1_frame) if exon_strand == '+' else c1_frame
            if c1_residue != 0:
                exon_seq_c1 = extract_seq(genome, exon_chrom, c1_cds['start'],c1_cds['end'], zero_based=0)[-c1_residue:]
            else:
                exon_seq_c1 = ""
            exon_seq = extract_seq(genome, exon_chrom, exon_s, exon_e, zero_based=1)
            exon_seq = exon_seq_c1 + exon_seq
            c2_frame = get_residue(c2_cds, c2_frame) if exon_strand == '-' else c2_frame
            if c2_frame != 0:
                exon_seq_c2 = extract_seq(genome, exon_chrom, c2_cds['start'],c2_cds['end'], 0)[:c2_frame]
            else:
                exon_seq_c2 = ""
            exon_seq = exon_seq + exon_seq_c2

            if exon_strand == '-':
                exon_seq = reverse_complement(exon_seq)
            
            # frame shift
            ptc_pos, ptc_seq = search_ptc(exon_seq)
            if ptc_pos != -1:
                impact_isoform.add('ORF-disrupting')
            
            impact_isoform.add(check_frame_shift(exon_seq))
            

           
        else:
            if c1_idx == -1 or c2_idx == -1:
                impact_isoform.add('False')
            else:
                impact_isoform.add('Not Found-{}_{}_{}_{}'.format(c1_idx,c2_idx, ts_s, ts_e))
    if len(impact_isoform) ==0:
        impact_isoform.add('ORF-preserving')
    impact_isoform = set(impact_isoform)
    if len(impact_isoform) == 1:
        if 'ORF-disrupting' in impact_isoform:
            return 'ORF-disrupting'
        else:
            return 'ORF-preserving'
    else:
        return 'ORF-disrupting-isoform'
    return ','.join(impact_isoform)


def get_dist_last_eej(ts_idx, gtf_ts, gtf_exon, gtf_ts_exon):

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
    gtf_ts = pd.read_csv(
        os.path.join(
            args.gtfdir, 'gencode.v31.primary_assembly.annotation.ts.tab'),
        sep='\t',
        dtype={'start': int, 'end': int, 'cds_start': int, 'cds_end': int, 'transcript_support_level':str})
    gtf_ts = gtf_ts[gtf_ts['gene_type'] == 'protein_coding']
    print("# of protein_coding transcripts: {}".format(gtf_ts.shape[0]))
    gtf_ts = gtf_ts[gtf_ts['transcript_support_level'].isin(['1',
                                                          '2'])]
    print("# of protein_coding tsl 1&2 transcripts: {}".format(gtf_ts.shape[0]))
     
    gtf_exon = pd.read_csv(
        os.path.join(
            args.gtfdir, 'gencode.v31.primary_assembly.annotation.exon.tab'),
        sep='\t', dtype={'start': int, 'end': int})
    gtf_exon = gtf_exon[gtf_exon['transcript_id'].isin(gtf_ts['transcript_id'])]
    print("# of exons: {}".format(gtf_exon.shape[0]))
    
    gtf_cds = pd.read_csv(
        os.path.join(
            args.gtfdir, 'gencode.v31.primary_assembly.annotation.cds.tab'),
        sep='\t', dtype={'start': int, 'end': int})
    gtf_cds = gtf_cds[gtf_cds['transcript_id'].isin(gtf_ts['transcript_id'])]
    print("# of cds: {}".format(gtf_cds.shape[0]))
    
    # cr_ts = build_cgranges(gtf_ts)
    cr_ts = build_cgranges(gtf_ts)
    
    # convert to dict
    gtf_ts = gtf_ts.to_dict('index')
    gtf_exon = gtf_exon.to_dict('index')
    gtf_cds = gtf_cds.to_dict('index')
    # build linking
    gtf_ts_exon = build_ts_exon_dict(gtf_exon)
    gtf_ts_cds = build_ts_exon_dict(gtf_cds)

    genome = twobitreader.TwoBitFile(args.genome)
    exon_target = pd.read_csv(
        args.bed, sep='\t',
        header=None,
        usecols=list(range(6)),
        names=['chrom', 'start', 'end', 'name', 'score', 'strand'])
    impact_flag_list = list()
    print(exon_target.head())
    for idx, exon in tqdm(exon_target.iterrows(), total=exon_target.shape[0]):
        # for idx, exon in exon_target.iterrows():
        
        impact_flag = track_ts(
            exon, genome=genome, cr_ts=cr_ts, gtf_ts=gtf_ts, gtf_exon=gtf_exon,
            gtf_ts_exon=gtf_ts_exon, gtf_ts_cds=gtf_ts_cds, gtf_cds=gtf_cds)
        impact_flag_list.append(impact_flag)
    exon_target['proteome_impact'] = impact_flag_list
    exon_target.to_csv(args.output, sep='\t', index=False)
