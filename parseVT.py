import os.path
import pandas as pd
import argparse
import swifter

"""
# Quality Score Column
# Score 1: Read coverage, based on actual reads
# Score 2: Read coverage, based on corrected reads
# Score 3: This score has been recicled to contain different information from release v2.2.2:
#   EX (except microexon module): raw total read counts supporting upstream inclusion, downstream inclusion and skipping (format INC1=INC2=EXC).
#   EX (microexon module): raw total read counts supporting inclusion and exclusion (format INC=EXC).
#   ALTD and ALTA: PSI-like value of the exon hosting the ALTD/ALTA event. This score is used to filter out events in compare based on the option --min_ALT_use.
#   IR (from v2.1.3): corrected number of intron body reads (in a sample of 200bp in the middle of the intron, or the whole intron if shorter), and the number of mappable position in that sample (maximum 151 positions) (format READS=POSITIONS).
#   Before v2.1.3: Read coverage, based on uncorrected reads mapping only to the reference C1A, AC2 or C1C2 splice junctions (similar values as per Score 1).
# Score 4: This score has different meaning depending on the type of AS event:
#   EX (except for microexon module): Imbalance of reads mapping to the inclusion splice junctions.
#       OK: the ratio between the total number of corrected reads supporting inclusion for splice junctions upstream and downstream the alternative exon is < 2.
#       B1: the ratio between the total number of corrected reads supporting inclusion for splice junctions upstream and downstream the alternative exon is > 2 but < 5.
#       B2: the ratio between the total number of corrected reads supporting inclusion for splice junctions upstream and downstream the alternative exon is > 5 (but none is 0).
#       B3: when the corrected reads for inclusion from one side is at least 15 and 0 for the other. Used to filter out events in compare when the option --noB3 is activated.
#       Bl/Bn: low (between 10 and 14)/no read coverage (between 1 and 9) for splice junctions supporting inclusion.
#   EX (microexon module): "na" (no information provided).
#   ALTD and ALTA: raw read counts for the specific splice site, for the all the splice sites of the event together (=total reads) and for those supporting skipping of the host exon. In versions early than v2.2.2: total reads for the event for all combinations or only for the reference acceptor (for ALTD) or donor (for ALTA).
#   IR: raw read counts mapping to the upstream exon-intron junction, downstream intron-exon junction, and exon-exon junction in the format EIJ=IEJ=EEJ). In versions earlier than v2.2.2, corrected counts were shown instead of raw read counts.
# Score 5:Complexity of the event
"""

def check_ce_thres(score):
    # Allowable coverage score (score 1) for CE
    cov_thres = ['SOK', 'OK', 'LOW']
    # Allowable balance score (score 4) for CE
    bal_thres = ['OK', 'B1']  
    score = score.split(',')
    check_res = score[0] in cov_thres and score[3] in bal_thres
    return check_res

def check_ir_thres(score):
    # Allowable coverage score (score 1) for IR
    cov_thres = 15
    bal_thres = 0.05
    score = score.split(',')
    nreads = sum([int(r) for r in score[3].split('=')])

    check_res = nreads >= cov_thres and float(score[4].split('@')[0]) > bal_thres
    return check_res

def check_other_thres(score):
    # Allowable coverage score (score 1) for "Alt5", "Alt3", "MIC"
    cov_thres = ['SOK', 'OK', 'LOW']
    score = score.split(',')
    check_res = score[0] in cov_thres 
    return check_res



def clean_AS(event, sample_list, min_prop=0.5):
    complex = event['COMPLEX']
    is_ce = complex in ['S', 'C1', 'C2', 'C3', 'ANN']
    is_ir = complex in ['IR']
    is_other = complex in ['Alt5', 'Alt3', 'MIC']
    check_complex = is_ce or is_ir or is_other
    if check_complex:
        if is_ce:
            check_list = event[sample_list].apply(check_ce_thres)
        if is_ir:
            check_list = event[sample_list].apply(check_ir_thres)
        if is_other:
            check_list = event[sample_list].apply(check_other_thres)
        check_prop = sum(check_list) / len(sample_list) >= min_prop
        return check_prop
    else:
        return check_complex



if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='''
                    Parse output from vast-tools
                    header of the default vast-tools align output:
                    col1: \t"GENE"\tCoressponding GENE
                    col2: \t"EVENT"\tEvent name
                    col3: \t"COORD"\tCoodinates
                    col4: \t"LENGTH"\tevent length
                    col5: \t"FullCO"\tFull coordinates
                    col6: \t"COMPLEX"\tType of even
                    col[7-:] \t"sample + score"
                    ''',
                    formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('-i', '--input', required=True, help='vast-tools align output')
    parser.add_argument('-o', '--output', help='filtered output')
    # parser.add_argument('-s', '--species',
    #                   choices=['Hg19', 'Hs2', 'Mm9', 'Mm10'])

    args = parser.parse_args()

    as_profile = pd.read_table(args.input)
    sample_list = as_profile.columns[7::2]
    print("Table containing {} sample(s) and {} event(s) for each sample:".format(len(sample_list), as_profile.shape[0]))
    print("including: {}".format(",".join(sample_list)))
    as_profile_filtered = as_profile[as_profile.swifter.apply(clean_AS, sample_list=sample_list,axis=1)]
    
    print("{} event(s) passed the filtering.".format(as_profile_filtered.shape[0]))
    if not args.output:
        file_name = os.path.splitext(args.input)
        args.output = file_name[0] + '_filtered' + file_name[1]
    as_profile_filtered.to_csv(args.output, sep="\t", index=False)