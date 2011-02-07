

def compareSampleLists(list1, list2):
    """ given two list of samplenames, compare element by element; return 1 if they are all the same 0 otherwise """
    if len(list1) != len(list2): return 0

    for i in range(0, len(list1) ):
        if list1[i] != list2[i]: return 0

    return 1

def get_vcfdataline(fh):
    """yield datalines """
    for line in fh:
        if '#' not in line: yield line.strip()
        

def yield_bedinterval(fh):
    """ yield tuple of (chr, start, end) zero-based half-open interval of vcf dataline """
    for line in fh:
        if '#' not in line:
            fields=split_vcfdataline( line.strip () )
            (chr,start)=fields[0:2]
            yield (chr, int(start)-1, int(start) )


def split_vcfdataline(line):
    """ spliti the fields of a vcf dataline and return in a list """
    return line.strip().split('\t')

def get_vcfdataline_passfilter(fh):
    """ yield list that has not been filtered """
    for line in fh:
        if '#' not in line:
            fields = line.strip().split('\t')
            if fields[6] == '.' or fields[6]=='PASS':
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
            if fields[6] == '.' or fields[6] == 'PASS':
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
    """return list of samples """
    for line in fh:
        if '#CHROM' in line:
            fields = line.strip().split('\t')
            samples = fields[9::]
            return samples

    
    

def stripGT(gt_string):
    """ strip a genotype string to only its GT field GT: => GT  """

    if gt_string == '.':
        return './.'
    else:
        gtformatfields=gt_string.split(':')
        return gtformatfields[0]


def returnAlleles(gt):
    """ return tuple (allele1, allele2) of stripped genotype GT field"""
    if ':' in gt:
        sys.stderr.write("strip string to only its GT field first!")
        return None
    if '/' in gt:
        (a1,a2)= gt.split('/')
        return (a1,a2)
    elif '|' in gt:
        (a1, a2) = gt.split('|')
        return (a1,a2)
    else:
        return None
        

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

def returnAlleles_phased(gt):
    """ return tuple of (allele1, allele2) of phased and stripped  GT string e.g. 0|0 returns (0,0) """
    if ':' in gt:
        sys.stderr.write("strip string to only its GT field first!")
        return None
    if '|' in gt:
        (a1,a2)= gt.split('|')
        return (a1,a2)
    else:
        return None


def compare_genotypes(g1, g2):
    """ given two lists of the form [ (sample, gt_string), ....  ], iterate thru and compare genotypes per individual    """
    """ return list of [ ( g1_genotype, g2_genotype, samplename) ] of missmatched genotypes"""

    compared_genotypes=[] # [ (g1_type, g2_type, sample) ... ]

    if len(g1) != len(g2):
        sys.stderr.write("cannot compare phased unphased genotypes; unequal size of lists!")
    else:
        for i in range( 0, len( g1 ) ):
            (p1,p2) = returnAlleles ( stripGT(g1[i][1] ) )
            (u1, u2) = returnAlleles ( stripGT(g2[i][1]) )

            g1_alleletype=typeofGenotype(p1,p2)
            g2_alleletype=typeofGenotype(u1,u2)
            #print g1[i], g1_alleletype
            #print g2[i], g2_alleletype
            compared_genotypes.append( (g1_alleletype, g2_alleletype, g1[i][0])  )
    return  compared_genotypes

def compare_phased_to_unphased(phased, unphased):
    """ give two lists of the form  [ (sample, gt_string), ....  ] in which one is unphased and theother phased iterate thru and compare the genotypes per individual  """

    unmatched_genotypes=[]
    if len(phased) != len(unphased):
        sys.stderr.write("cannot compare phased unphased genotypes; unequal size of lists!")
    else:
        for i in range( 0, len( phased ) ):
            if stripGT(unphased[i][1]) != '.':
                (p1,p2) = returnAlleles_phased ( stripGT(phased[i][1] ) )
                (u1, u2) = returnAlleles_unphased ( stripGT(unphased[i][1]) )
                if getNonRefDosage(p1,p2) != getNonRefDosage(u1,u2):
                    unmatched_genotypes.append( (phased[i][1], unphased[i][1], unphased[i][0] )    )
    return  unmatched_genotypes

    
def getNonRefDosage(allele1,allele2):
    """ return the number of non-ref alleles given the two alleles from a genotype """

    dosage=0
    if '1' in allele1: dosage+=1
    if '1' in allele2: dosage+=1

    return dosage

def typeofGenotype(allele1, allele2):
    """ I really should be a python version of a typedef here, but dont know how 
        hom_ref =1 het =2 hom_nonref=3 no_call=4                              """

    if allele1 == '0' and allele2 == '0': return 0

    if  allele1 == '0' and allele2== '1': return 1
    if allele1 =='1' and allele2 == '0': return 1

    if allele1== '1' and allele2== '1': return 2

    if allele1== '.' or allele2 == '.': return 3
    if allele1 == '.' and allele2 == '.': return 3