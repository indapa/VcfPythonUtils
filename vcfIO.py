def get_vcfdataline(fh):
    """yield datalines """
    for line in fh:
        if '#' not in line: yield line.strip()
        


def get_vcfdataline_passfilter(fh):
    """ yield list that has not been filtered """
    for line in fh:
        if '#' not in line:
            fields = line.strip().split('\t')
            if fields[6] == '.':
                yield line.strip().split('\t')

def get_vcfdatafields(fh):
    """ yield data fields """
    for line in fh:
        if '#' not in line:
            fields = line.strip().split('\t')
            yield fields[0:9]

def get_vcfgenotypes(fh):
    """ yield genotypes from vcf dataline """
    for line in fh:
        if '#' not in line:
            fields = line.strip().split('\t')
            yield fields[9::]


def get_vcftuples_keep(fh, keepList):
    """yield  (chrom, pos, ref, alt, [ (sample,gt), ...., (sample,gt) ] ) but only for those samples in the keepList"""
    samples = get_vcfsamples(fh)
    for line in fh:
        if '#' not in line:
            fields = line.strip().split('\t')
            tuple = zip(samples,fields[9::])
            filtered_tuple= [x for x in tuple if x[0] in keepList]
            yield(fields[0], fields[1], fields[3], fields[4], filtered_tuple)

def get_vcftuple_passfilter(fh):
    """yield tuple with (chrom, pos, ref, alt, [ (sample,gt), ...., (sample,gt) ] )  of datalines that were not filtered"""
    samples = get_vcfsamples(fh)
    for line in fh:
        if '#' not in line:
            fields = line.strip().split('\t')
            if fields[6] == '.':
                t=zip(samples,fields[9::])
                yield (fields[0],fields[1],fields[3], fields[4], t)



def get_vcftuples(fh):
    """yield tuple with (chrom, pos, ref, alt, [ (sample,gt), ..., (sample,gt) ] ) of datalines regarless of filter field"""
    samples = get_vcfsamples(fh)
    for line in fh:
        if '#' not in line:
            fields = line.strip().split('\t')
            t=zip(samples,fields[9::])
            yield (fields[0], fields[1], fields[3], fields[4], t)


def get_vcfsamples_keep(fh, keepList):
    """ yield a list of samples that are in the keepList """
    for line in fh:
        line.strip()
        if '#CHROM' in line:
            fields = line.strip().split('\t')
            samples = fields[9::]
            kept=[s for s in samples if s in keepList]
            return kept

def get_vcfsamples(fh):
    """yield list of samples """
    for line in fh:
        line.strip()
        if '#CHROM' in line:
            fields = line.strip().split('\t')
            samples = fields[9::]
            return samples

def stripGT(gt_string):
    """ strip a genotype string to only its GT field GT: => GT  """
    gtformatfields=gt_string.split(':')
    return gtformatfields[0]

def returnAlleles_unphased(gt):
    """ return tuple of (allele1, allele2) of unphased and stripped  GT string e.g. 0/0 returns (0,0) """
    if ':' in gt:
        sys.stderr.write("strip string to only its GT field first!")
        return None
    if '/' in gt:
        (a1,a2)= gt.split('/')
        return (a1,a2)
    else:
        return None


