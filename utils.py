import pandas
import swifter
import re
def reverse_complement(seq):
    try:
        assert set(seq).issubset(set("ACGTN"))
    except AssertionError:
        raise ValueError(
            "Sequence {} contains invalid values: {}"
            .format(seq, set(seq) - set('ACGTN'))
        )
    return ''.join('TGCAN'['ACGTN'.index(s)] for s in seq.upper()[::-1])

def parse_gtf_attrs(feature):
    feature_type = feature['feature']
    attrs = feature['attribute']
    regex = r"(\S+\s\S+;)"
    attrs_dict = dict()
    for attr in re.finditer(regex,attrs):
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
    if feature_type == 'transcript':
        feature['gene_id'] = attrs_dict['gene_id']
        feature['gene_type'] = attrs_dict['gene_type']
        feature['gene_name'] = attrs_dict['gene_name']
        feature['transcript_id'] = attrs_dict['transcript_id']
        if 'transcript_support_level' in attrs_dict:
            feature['transcript_support_level'] = attrs_dict['transcript_support_level']
        else:
            feature['transcript_support_level'] = 'NA'
    if feature_type == 'exon':
        feature['gene_id'] = attrs_dict['gene_id']
        feature['gene_type'] = attrs_dict['gene_type']
        feature['gene_name'] = attrs_dict['gene_name']
        feature['transcript_id'] = attrs_dict['transcript_id']
        feature['exon_id'] = attrs_dict['exon_id']
        feature['exon_number'] = int(attrs_dict['exon_number'])
    if feature_type == 'CDS':
        feature['gene_id'] = attrs_dict['gene_id']
        feature['gene_type'] = attrs_dict['gene_type']
        feature['gene_name'] = attrs_dict['gene_name']
        feature['transcript_id'] = attrs_dict['transcript_id']
        feature['exon_id'] = attrs_dict['exon_id']
        feature['exon_number'] = int(attrs_dict['exon_number'])

    return feature



def parse_gtf(gtf_file):
    gtf = pandas.read_csv(gtf_file, sep='\t', comment='#', names=['seqname','source','feature','start','end','score','strand','frame','attribute'])
    gtf_gene = gtf[gtf['feature']=='gene'].swifter.apply(parse_gtf_attrs,axis=1)
    gtf_transcript = gtf[gtf['feature']=='transcript'].swifter.apply(parse_gtf_attrs,axis=1)
    gtf_exon = gtf[gtf['feature']=='exon'].swifter.apply(parse_gtf_attrs,axis=1)
    gtf_CDS = gtf[gtf['feature']=='CDS'].swifter.apply(parse_gtf_attrs,axis=1)
    return gtf_gene, gtf_transcript, gtf_exon,gtf_CDS
