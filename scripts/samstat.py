import argparse

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Extract samtools stat information')
    parser.add_argument('-i', '--input', help='input file')
    args = parser.parse_args()

    with open(args.input) as stat_file:
        count_sum = 0
        percentile_sum = 0
        for line in stat_file:
            if line.startswith('GCF') or line.startswith('GCL'):
                tag, percentile, count = line.rstrip().split()
                percentile = float(percentile)
                count = int(count)
                count_sum = count_sum + count
                percentile_sum = percentile_sum + count * percentile
        print("{}\t{}\t{}".format(args.input, percentile_sum, count_sum))
