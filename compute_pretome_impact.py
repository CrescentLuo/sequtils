import argparse
import pybedtools
import bisect
import twobitreader
from operator import itemgetter, attrgetter

def get_args():
    parser = argparse.ArgumentParser(description='Compute the impact of novel exons')
    parser.add_argument('-g', '--gtf', required=True, help='input annotation file in GTF format')
    parser.add_argument('-b', '--bed', required=True, help='table of novel exons in BED format')
    parser.add_argument('--genome', required=True, help='genome sequence file in 2bit format')
    args = parser.parse_args()
    return args

class AnnoGTF:
    def __init__(self):
        self.gene_set = set()
        self.gene_id_dict = dict()
        self.gene_ts_dict = dict()
        self.transcript_exon_dict = dict()
        self.transcript_cds_dict = dict()
        self.gene_set = dict()
        self.transcript_set = dict()
        self.exon_set = dict()

    def init(self,file):
        iterms = pybedtools.BedTool(file)
        gene_set = iterms.filter(lambda iterm: iterm.fields[2] == 'gene')
        transcript_set = iterms.filter(lambda iterm: iterm.fields[2] == 'transcript')
        cds_set = iterms.filter(lambda iterm: iterm.fields[2] == 'CDS')
        exon_set = iterms.filter(lambda iterm: iterm.fields[2] == 'exon')
        

        for gene in gene_set:
            gene_id = gene.attrs['gene_id']
            gene_symbol = gene.attrs['gene_name']
            self.gene_id_dict[gene_id] = gene_symbol
            self.gene_id_dict[gene_symbol] = gene_id
            self.gene_ts_dict[gene_id] = set()
            self.gene_set[gene_id] = gene

        for ts in transcript_set:
            ts_id = ts.attrs['transcript_id']
            gene_id = ts.attrs['gene_id']
            self.gene_ts_dict[gene_id].add(ts_id)
            self.transcript_exon_dict[ts_id] = set()
            self.transcript_set[ts_id] = ts
            self.transcript_cds_dict[ts_id] = list()

        for exon in exon_set:
            ts_id = exon.attrs['transcript_id']
            exon_id = exon.attrs['exon_id']
            self.transcript_exon_dict[ts_id].add(exon_id)
            self.exon_set[exon_id] = exon
        
        for cds in cds_set:
            ts_id = cds.attrs['transcript_id']
            self.transcript_cds_dict[ts_id].append(cds)
            
        
        return self
    
    def get_ts(self, gene_id):
        try:
            return self.gene_ts_dict[gene_id]
        except:
            return
    

def read_gtf(gtf_file):
    gtf = AnnoGTF().init(file=gtf_file)
    #subset = a.filter(lambda b: b.chrom == 'chr1' and b.start < 150)
    
    return gtf

def read_exon(exon_file):
    exon = pybedtools.BedTool(exon_file)
    return exon

def search_stop_codon(seq, lo=0):
    # stop codon TAG TAA TGA
    stop_codon = set("TAG", "TAA", "TGA")
    for i in range(0,len(seq),3):
        codon = seq[i:i+2]
        if codon in stop_codon:
            return i, codon
    return -1,"NA"




def track_ts(exon, gtf, genome):
    gene = exon.name.split("_")[1]
    print(gene)
    gene_id = gtf.gene_id_dict[gene]
    ts_set = gtf.get_ts(gene_id)
    if exon.strand == '+':
        exon_seq = genome[exon.chrom][exon.start:exon.end]
    else:

    for ts_id in ts_set:
        ts = gtf.transcript_set[ts_id]
        cds_set_sort = list()
        if exon.start >= ts.start and exon.end <= ts.end:
            cds_set = gtf.transcript_cds_dict[ts_id]
            cds_set = sorted(cds_set)
            pos = -1
            cds_seq = ''
            for idx, cds in enumerate(cds_set):
                if cds > exon:
                    pos = idx
                    residue = len(cds_seq) % 3
                    if residue:
                        residue_seq = cds_seq[-residue::]
                        exon_seq_with_overhang = residue_seq + exon_seq[::residue-3]
                        
                    else:
                    
                    pos, codon = search_stop_codon(exon_seq_with_overhang)
                    if pos != -1:
                    # disruption
                        print pos, codon
                        return
                     


                else:
                    cds_seq = cds_seq + genome[cds.chrom][cds.start:cds.end]
                    

        
        


if __name__ == '__main__':
    args = get_args()
    gtf = read_gtf(args.gtf)
    exons = read_exon(args.bed)
    genome = twobitreader.TwoBitFile(args.genome)
    for e in exons:
        track_ts(e, gtf, genome)
    
    