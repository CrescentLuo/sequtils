import argparse
import swifter
import pandas as pd
import pysam
from utils import reverse_complement

def seq_gc_content(seq):
  return sum(s in {'G', 'C'} for s in seq) / len(seq)

def calculate_gc_percentage(interval, genome,extend=0):
    chrom = interval['chrom']
    start = int(interval['start'])
    end = int(interval['end'])
    strand = interval['strand']
    genome = pysam.FastaFile(genome)
    if strand == '+' or strand == '.':
        mid_region_seq = genome.fetch(chrom,start,end).upper()
        interval['GC_exon'] = seq_gc_content(mid_region_seq)
        if extend != 0:
            up_intron_seq = genome.fetch(chrom,start-extend,start).upper()
            dn_intron_seq = genome.fetch(chrom,end, end+extend).upper()
            interval['GC_up_intron'] = seq_gc_content(up_intron_seq)
            interval['GC_dn_intron'] = seq_gc_content(dn_intron_seq)
    else:
        mid_region_seq = reverse_complement(genome.fetch(chrom,start,end).upper())
        interval['GC_exon'] = seq_gc_content(mid_region_seq)
        if extend != 0:
            up_intron_seq = reverse_complement(genome.fetch(chrom,end, end+extend).upper())
            dn_intron_seq = reverse_complement(genome.fetch(chrom,start-extend,start).upper())
            interval['GC_up_intron'] = seq_gc_content(up_intron_seq)
            interval['GC_dn_intron'] = seq_gc_content(dn_intron_seq)
    return interval

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='compute GC content of interested intervals')
    parser.add_argument('-b', '--bed', required=True, help='Input file in BED format')
    parser.add_argument('-g', '--genome', required=True, help='Reference genome file')
    parser.add_argument('-e', '--extend', type=int, default=0)
    parser.add_argument('-o', '--output', default='./GC_content_output.tab',
                        help='Output GC content tab file')
    args = parser.parse_args()
    gc_record = pd.read_csv(args.bed, sep='\t', names=['chrom','start','end','name','score','strand'])
    gc_record = gc_record.swifter.apply(calculate_gc_percentage, genome=args.genome, extend=args.extend,axis=1)
    gc_record.to_csv(args.output, sep='\t', index=False)