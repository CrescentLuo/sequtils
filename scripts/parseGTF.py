import pandas as pd
import re
import pysam
from utils import reverse_complement
import pathlib

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
    elif feature_type == 'transcript':
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
    elif feature_type == 'exon':
        feature['gene_id'] = attrs_dict['gene_id']
        feature['gene_type'] = attrs_dict['gene_type']
        feature['gene_name'] = attrs_dict['gene_name']
        feature['transcript_id'] = attrs_dict['transcript_id']
        feature['exon_id'] = attrs_dict['exon_id']
        feature['exon_number'] = int(attrs_dict['exon_number'])
    elif feature_type == 'CDS':
        feature['gene_id'] = attrs_dict['gene_id']
        feature['gene_type'] = attrs_dict['gene_type']
        feature['gene_name'] = attrs_dict['gene_name']
        feature['transcript_id'] = attrs_dict['transcript_id']
        feature['exon_id'] = attrs_dict['exon_id']
        feature['exon_number'] = int(attrs_dict['exon_number'])

    return feature



def parse_gtf(gtf_file):
    gtf = pd.read_csv(
        gtf_file, sep='\t', comment='#',
        names=['seqname', 'source', 'feature', 'start', 'end', 'score',
                'strand', 'frame', 'attribute'])
    gtf_gene = gtf[gtf['feature'] == 'gene'].swifter.apply(
        parse_gtf_attrs, axis=1)
    gtf_ts = gtf[gtf['feature'] == 'transcript'].swifter.apply(parse_gtf_attrs, axis=1)
    gtf_exon = gtf[gtf['feature'] == 'exon'].swifter.apply(parse_gtf_attrs, axis=1)
    gtf_cds = gtf[gtf['feature'] == 'CDS'].swifter.apply(parse_gtf_attrs, axis=1)
    
    gtf_ts = assign_cds_to_ts(gtf_ts, gtf_cds)

    return gtf_gene, gtf_ts, gtf_exon, gtf_cds

    
def concat_interval_seq(ts_id, ts_dict, interval_set,genome_seq):
    concat_seq = ''
    strand = '+'
    #print(ts_dict[ts_id])
    for interval_id in ts_dict[ts_id]:
        chrom = interval_set[interval_id]['seqname']
        if interval_set[interval_id]['strand'] == '-':
            strand = '-'
        start = int(interval_set[interval_id]['start'])
        start = start - 1
        end = int(interval_set[interval_id]['end'])
        concat_seq = concat_seq + genome_seq.fetch(chrom,start,end)
    if strand == '-':
        return reverse_complement(concat_seq)
    return concat_seq

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

def build_ts_dict(gtf_cds):
    ts_dict = dict()
    for idx in gtf_cds:
        transcript_id = gtf_cds[idx]['transcript_id']
        if transcript_id not in ts_dict:
            ts_dict[transcript_id] = list()
        ts_dict[transcript_id].append(idx)
    
    return ts_dict

def extract_annotation(annotation_dir):
    gtf_gene_annot = sorted(pathlib.Path(annotation_dir).glob("*.gene.tab"))[0]
    gtf_gene = pd.read_csv(gtf_gene_annot, sep='\t')
    gtf_ts_annot = sorted(pathlib.Path(annotation_dir).glob("*.ts.tab"))[0]
    gtf_ts = pd.read_csv(gtf_gene_annot, sep='\t')
    return gtf_gene, gtf_ts