import argparse
import pybedtools

def get_flanking_region(feature, upstream, width):
    if upstream:
        if feature.strand == '+':
            feature.end = feature.start
            feature.start = feature.start - width
        if feature.strand == '-':
            feature.start = feature.end
            feature.end = feature.end + width
    else:
        if feature.strand == '+':
            feature.start = feature.end
            feature.end = feature.end + width
        if feature.strand == '-':
            feature.end = feature.start
            feature.start = feature.start - width
    return feature


if __name__ == '__main__':
    argparser = argparse.ArgumentParser(
        description="extract flanking regions of target regions from input BED file")
    argparser.add_argument('-i', '--input', required=True, help='input BED file')
    argparser.add_argument('-w', '--width', type=int, default=200, help='width of flanking region')
    argparser.add_argument('-d', '--direct', choices=['up', 'dn', 'both'], help="direction of target region")
    argparser.add_argument('-p', '--prefix', default="out", help="output file")
    args = argparser.parse_args()

    input_bed = pybedtools.BedTool(args.input)
    if args.direct == 'up':
        flank_region_up = input_bed.each(get_flanking_region, upstream=True, width=args.width)
        flank_region_up.saveas("{}_up{}.bed".format(args.prefix, args.width))
    elif args.direct == 'dn':
        flank_region_dn = input_bed.each(get_flanking_region, upstream=False, width=args.width)
        flank_region_dn.saveas("{}_dn{}.bed".format(args.prefix, args.width))
    else:
        flank_region_up = input_bed.each(get_flanking_region, upstream=True, width=args.width)
        flank_region_dn = input_bed.each(get_flanking_region, upstream=False, width=args.width)
        flank_region_up.saveas("{}_up{}.bed".format(args.prefix, args.width))
        flank_region_dn.saveas("{}_dn{}.bed".format(args.prefix, args.width))
