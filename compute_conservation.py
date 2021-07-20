import argparse
import pyBigWig as pyBW
import numpy as np
import pandas as pd
import swifter

# Phastcons BigWig at : 
# http://hgdownload.cse.ucsc.edu/goldenPath/hg38/phastCons7way/hg38.phastCons7way.bw
# http://hgdownload.cse.ucsc.edu/goldenPath/hg38/phyloP7way/hg38.phyloP7way.bw

def arg_parser():
    
    return args

def get_phastcons(interval, phastcons_bw, nBins=10):

    chrom = interval['chrom']
    start = int(interval['start'])
    end = int(interval['end'])
    interval['phastcons_bin'+'(bins={})'.format(nBins)] = ','.join([str(s) for s in phastcons_bw.stats(chrom,start,end,nBins=nBins)])
    interval['phastcons_mean'] = phastcons_bw.stats(chrom,start,end)[0]
    return interval


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Calculation of PhastCons of genes using bigwig file")
    parser.add_argument('-p', '--phastcons', required=True, help='Phast conservation in bigwig format')
    parser.add_argument('-b', '--bed', required=True, help='Interested gene genomic location in BED format')
    parser.add_argument('-n', '--nBins', type=int, default=10, help='number of bins for phastcons calculation')
    parser.add_argument('-o', '--output', default='./phastcons_output.tab',
                        help='Output MaxEntScan score tab file')
    args = parser.parse_args()

    phastcons_record = pd.read_csv(args.bed, sep='\t', names=['chrom','start','end','name','score','strand'])
    with pyBW.open(args.phastcons) as phastcons_bw:
        phastcons_record = phastcons_record.swifter.apply(get_phastcons, phastcons_bw=phastcons_bw, axis=1)
        phastcons_record.to_csv(args.output, sep='\t', index=False)
        