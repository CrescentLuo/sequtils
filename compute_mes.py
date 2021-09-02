import argparse
from os import replace
import maxentpy
import pysam
from utils import reverse_complement
from maxentpy.maxent import load_matrix5, load_matrix3
from maxentpy import maxent_fast
import swifter

import pandas as pd

def parse_args():
    parser = argparse.ArgumentParser(
        description='MaxEntScan score for bed files')
    parser.add_argument('-b', '--bed', required=True,
                        help='Input interval file in bed format')
    parser.add_argument('-g', '--genome', required=True,
                        help='Genome fasta file(Indexed)')
    parser.add_argument('-o', '--output', default='./mes_output.tab',
                        help='Output MaxEntScan score tab file')
    args = parser.parse_args()
    return args

def compute_mes(interval, matrix5, matrix3, genome):
    genome = pysam.FastaFile(genome)
    # 5ss 3bases in exon and 6 bases in intron
    # 3ss 20 bases in intron and 3 bases in exon
    if interval['strand'] == '+':
        seq5 = genome.fetch(interval.chrom, interval.end-3, interval.end + 6).upper()
        seq3 = genome.fetch(interval.chrom, interval.start - 20, interval.start + 3).upper()
    else:
        seq5 = reverse_complement(genome.fetch(
                interval.chrom, interval.start - 6, interval.start + 3).upper())
        seq3 = reverse_complement(genome.fetch(
                interval.chrom, interval.end - 3, interval.end + 20).upper())
    name_format_str = '{seq5}:{mes5}|{seq3}:{mes3}'
    if set(seq5).issubset('ACGT') and set(seq3).issubset('ACGT'):
            mes5 = maxent_fast.score5(seq5, matrix=matrix5)
            mes3 = maxent_fast.score3(seq3, matrix=matrix3)
            interval['seq5'] = seq5
            interval['mes5'] = mes5
            interval['seq3'] = seq3
            interval['mes3'] = mes3
            
    else:
        interval['seq5'] = seq5
        interval['mes5'] = 'NA'
        interval['seq3'] = seq3
        interval['mes3'] = 'NA'
    return interval

if __name__ == '__main__':
    args = parse_args()
    mes_record = pd.read_csv(args.bed, sep='\t', names=['chrom','start','end','name','score','strand'])

    # load maxent mamtrix
    matrix5 = load_matrix5()
    matrix3 = load_matrix3()

    mes_record = mes_record.swifter.apply(compute_mes, matrix5=matrix5, matrix3=matrix3, genome=args.genome, axis=1)
    mes_record['mes5'] = pd.to_numeric(mes_record['mes5'])
    mes_record['mes3'] = pd.to_numeric(mes_record['mes3'])
    mes_mean = mes_record[['mes5', 'mes3']].mean()
    mes_cov = mes_record[['mes5', 'mes3']].cov()
    print(mes_mean)
    print(mes_cov)
    mes_record.to_csv(args.output, sep='\t', index=False)
