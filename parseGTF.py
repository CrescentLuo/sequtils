
class GtfFeature():
    sequence = None
    chrom = None
    source = None
    feature_type = None
    start = None
    end = None
    score = '.'
    strand = None
    phase = None
    attributes = dict()

    def __init__(self, feature):
        if isinstance(feature, list):
            self.sequence = feature[0]
            self.chrom = slef.sequence
            self.source = feature[1]
            self.feature_type = feature[2]
            self.start = int(feature[3])
            self.end = int(feature[4])
            self.score = feature[5]
            self.strand = feature[6]
            self.phase = feature[7]
            self.attributes = dict(attr.lstrip().replace('"', "").split(
                ' ') for attr in sline[8].rstrip(';').split(';'))


class GtfParser():
    gtf_file = None
    gtf_instance = None
    iterms = list()
    gene_to_transcripts = dict()
    transcript_to_exons = dict()
    

    def __init__(self,
                 fn=None,
                 ):
        # init gtf annotation
        self.gtf_file = fn
        print(self.gtf_file)
        self.i
        with open(self.gtf_file) as gtf:
            for line in gtf:
                if line.startswith('#'):
                    continue
                feature = line.rstrip.split('\t')
                GtfFeature(feature)

        self.iterms = [i for i in self.gtf_instance]

    def head(self, n=10):
        for i, interval in enumerate(self.iterms):
            if i == n:
                break
            print(interval)
            # print("\t".join(interval))

    def filter(self, func, *args, **kwargs):
        """
        Filter features by user-defined function.

        Takes a function *func* that is called for each feature in the
        `BedTool` object and returns only those for which the function returns
        True.

        *args and **kwargs are passed directly to *func*.

        Returns a streaming BedTool; if you want the filename then use the
        .saveas() method.

        #>>> a = pybedtools.example_bedtool('a.bed')
        #>>> subset = a.filter(lambda b: b.chrom == 'chr1' and b.start < 150)
        #>>> len(a), len(subset)
        (4, 2)

        so it has extracted 2 records from the original 4.

        """
        # TODO: write the return value of filter
        
        # return BedTool((f for f in self if func(f, *args, **kwargs))
