import sys
import re
import math

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
                yield line.strip()

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
    
def getFormatfield(gt_string, index):
    """ return the ith field of a genotype format string specified by index"""
    gtformatfields=gt_string.split(':')
    return gtformatfields[index]

def stripGT(gt_string):
    """ strip a genotype string to only its GT field GT: => GT  """

    if gt_string == '.':
        return './.'
    else:
        gtformatfields=gt_string.split(':')
        return gtformatfields[0]

def yieldGQ (g1, formatstr):
    """ given a list of [ (sample, gt_string) ] parse out the GQ and yield the genotype quality for the genotype """
    formatfields=formatstr.split(':')
    if 'GQ' in formatfields:
        gq_index = formatfields.index('GQ')
    else:
        sys.stderr.write("GQ not in genotype format string, cannot return genotype quality!\n")
        exit(1)
    for (sample, gtstring) in g1:
        (allele1, allele2) = returnAlleles( stripGT(gtstring)  )
        if (allele1 == '.' or allele2 == '.'): continue
        gq=getFormatfield(gtstring, gq_index)
        yield gq
    


def returnAlleles(gt):
    """ return tuple (allele1, allele2) of stripped genotype GT field"""

    if ':' in gt:
        sys.stderr.write("strip string to only its GT field first!")
        return None
    elif '/' in gt:
        (a1,a2)= gt.split('/')
        return (a1,a2)
    elif '|' in gt:
        (a1, a2) = gt.split('|')
        return (a1,a2)
    elif '.' in gt:
        return ('.', '.')
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

def getCalledGenotypeCount (g):
    """ return the numebr of called genotpes given a list  [ (sample, gt_string) ... ] return the number of called genotypes """
    total=0
    for (sample, genotype) in g:
        (p1,p2) = returnAlleles ( stripGT( genotype ) )

        if typeofGenotype(p1,p2) != 3:
            total+=1

    return total


def calMaf (g):
    """ given a list in the form [ (sample, gt_string) ] determine the minor allele freq= alt/2*called genotypes """
    total=0
    alt_count=0
    for (sample, genotype) in g:
        (p1,p2) = returnAlleles ( stripGT( genotype ) )
        
        alleletype=typeofGenotype(p1,p2)

        if alleletype == 3:
            continue
        total+=1
        if alleletype ==1 : 
            alt_count+=1
        elif alleletype == 2:
            alt_count +=2
        else:
            pass
    maf = float(alt_count) / float(2*total)
    #print float(2*total)
    #print p1, p2, sample, maf, alt_count, total
    return  ( round(maf,3) )


def return_expected_genotype_counts ( g1 ):
    """ given a list in the form [ (sample, gt_string), ... ] return expected counts of 3 genotypes """

    theta = calMaf(g1)

    pAA =  float(1-theta) * float(1-theta) * float( getCalledGenotypeCount(g1) )
    pAB =  float( 2*theta*(1-theta) ) * float ( getCalledGenotypeCount(g1) )
    pBB =  float( theta) * float (theta) * float ( getCalledGenotypeCount(g1) )

    return (round(pAA,3), round(pAB,3), round(pBB,3) )


def return_observed_genotype_counts( g1 ):
    """ given a list in the form [ (sample, gt_string), ... ] return the counts of the 3 types of genotypes """
    homoz_ref=0
    het=0
    homoz_nonref=0
    for i in range (0, len( g1 ) ):
        (p1,p2) = returnAlleles ( stripGT( g1[i][1] ) )
        if typeofGenotype(p1,p2) == 3:
            continue
        else:
            alleletype=typeofGenotype(p1,p2)
            #print g1[i][0], alleletype, g1[i][1]
            if alleletype==0:
                homoz_ref+=1
            if alleletype==1:
                het+=1
            if alleletype==2:
                homoz_nonref+=1
    return( homoz_ref,  het, homoz_nonref )

def hweLRT ( g1 ):
    """ return the LRT for HWE """

    #2 * sum (obs * log(obs/exp)) = LRT

    lrt=0

    obs = (pAAobs, pABobs, pBBobs ) = return_observed_genotype_counts ( g1 )
    exp = (pAAexp, pABexp, pBBexp) =  return_expected_genotype_counts ( g1 )

    for i in range(0, len(obs)):
        if obs[i] == 0:
            lrt+=0
        else:
            lrt+= float(obs[i]) * math.log(obs[i]/exp[i])

    lrt *= 2
    return  lrt


def gq_calibration ( g1, truth, formatstr):
    """ given two lists in the form [ (sample, gt_string), ... ] and the second list are 'truth' genotypes return a list of [ (gq, 0|1) ] were 1
    indicates the genotype in g1 matched the truth and 0 indicates it did not. Note, 'GQ' needs to be in format string inorder to calibrate """

    formatfields = formatstr.split(':')
    if 'GQ' in formatstr:
        gq_index = formatfields.index('GQ')
        
    else:
        sys.stderr.write("genotype format string doesn't contain GQ, cannot perform genotype quality calibration!\n")
        exit(1)

    gq_calibration = []
    
    if len(g1) != len(truth):
        sys.stderr.write("gq_calibration: cannot compare  genotypes; unequal size of lists!\n")
        exit(1)
    else:
        for i in range(0, len(g1) ):
            if g1[i][0] != truth[i][0]:
                sys.stderr.write("gq_calibration:samples don't match in comparison\n")
                exit(1)
            (p1,p2) = returnAlleles ( stripGT( g1[i][1] ) )
            if typeofGenotype(p1,p2) == 3: continue
            
            gq=getFormatfield(g1[i][1], gq_index)
            
            (u1, u2) = returnAlleles ( stripGT( truth[i][1]) )
            if typeofGenotype(u1,u2) == 3: continue

            correctly_genotyped=0

            g1_alleletype=typeofGenotype(p1,p2)
            g2_alleletype=typeofGenotype(u1,u2)

            if g1_alleletype != 3 and g2_alleletype !=3 and ( g1_alleletype != g2_alleletype ):
                correctly_genotyped=0
            elif g1_alleletype != 3 and g2_alleletype !=3 and ( g1_alleletype == g2_alleletype ):
                correctly_genotyped=1
                #print g1_maxprob, g2[i][1], correctly_imputed
            else:
                pass

            gq_calibration.append( (gq, correctly_genotyped) )
            
    return gq_calibration

def posterior_imputed_gprob_calibration ( g1, g2, formatstr):
    """ given two lists of the form [ (sample, gt_string), ... ] and the second one is the 'truth' genotypes
        return a list [ (max_gprob, 0|1), ... ] where 1 is the imputed genotype matched the truth 0 is it did not
        note the list returned will only be as long as how many imputed gentypes there were for that marker  """

    formatfields=formatstr.split(':')
    if 'GPROB' in formatstr and 'OG' in formatstr:
        gprobs_index=formatfields.index('GPROB')
        og_index=formatfields.index('OG')
    else:
        sys.stderr.write("genotype format doesn't contain GPROB and/or OG , cannot perform posterior genotype calibration!\n")
        exit(1)
    imputed_genotype_calibration=[]

    if len(g1) != len(g2):
        sys.stderr.write("posterior_imputed_gprob_calibration: cannot compare  genotypes; unequal size of lists!\n")
    else:
        for i in range( 0, len( g1 ) ):

            if g1[i][0] != g2[i][0]:
                sys.stderr.write("samples don't match in genotypes comparison!")
                exit(1)

            #get the original genotype and make sure it was originally missing data (only want to compare imputed genotypes)
            (og1, og2) = returnAlleles ( getFormatfield( g1[i][1], og_index ) )
            if (typeofGenotype(og1,og2) != 3):
                continue

            #get the genotyep proba of the imputed genotype
            g1_maxprob= max ( getFormatfield( g1[i][1], gprobs_index ).split(';') )
            
            
            (p1,p2) = returnAlleles ( stripGT( g1[i][1] ) )
            (u1, u2) = returnAlleles ( stripGT( g2[i][1]) )

            

            g1_alleletype=typeofGenotype(p1,p2)
            g2_alleletype=typeofGenotype(u1,u2)
            correctly_imputed=0

            if g1_alleletype != 3 and g2_alleletype !=3 and ( g1_alleletype != g2_alleletype ):
                correctly_imputed=0
                #print  g1_maxprob, g1[i][1], g2[i][1], correctly_imputed
                
            elif g1_alleletype != 3 and g2_alleletype !=3 and ( g1_alleletype == g2_alleletype ):
                correctly_imputed=1
                #print g1_maxprob, g2[i][1], correctly_imputed
            else:
                pass

            imputed_genotype_calibration.append( (g1_maxprob, correctly_imputed) )
            

    return imputed_genotype_calibration


def compare_nonimputed_genotypes( g1, g2, formatstr, thresh):
    """ given two lists of the form [ (sample, gt_string), ... ] iterate and compare *non-imputed* genotypes in g1 to genotypes in g2  """
    """ to find out which genotypes in g1 are imputed the FORMAT string must contain the OG tag and the OG must be a nocall """
    """ return list of [ ( g1_alleletype, g2_alleletype, samplename) ] where alleletype is [0, homref; 1, het; 2 hom_nonref; 3, nocall  """

    formatfields=formatstr.split(':')
    if 'GPROB' in formatstr and 'OG' in formatstr:
        gprobs_index=formatfields.index('GPROB')
        og_index=formatfields.index('OG')
    else:
        sys.stderr.write("genotype format doesn't contain GPROB/OG, cannot determine if genotype was imputed!")
        exit(1)
    compared_imputed_genotypes=[]

    if len(g1) != len(g2):
        sys.stderr.write("cannot compare phased unphased genotypes; unequal size of lists!")
    else:
        for i in range( 0, len( g1 ) ):

            if g1[i][0] != g2[i][0]:
                sys.stderr.write("samples don't match in genotypes comparison!")
                exit(1)

            #check if the original genotype in g1 was imputed; if not pass on the comparison!
            (og1, og2) = returnAlleles ( getFormatfield( g1[i][1], og_index ) )
            if (typeofGenotype(og1,og2) == 3):
                continue
            g1_maxprob= max ( getFormatfield( g1[i][1], gprobs_index ).split(';') )
            g1_maxprob=float(g1_maxprob)
            #posterior prob of imputed genotype did not meet threshold!
            if g1_maxprob <=float(thresh):
                continue
            (p1,p2) = returnAlleles ( stripGT( g1[i][1] ) )
            (u1, u2) = returnAlleles ( stripGT( g2[i][1]) )

            g1_alleletype=typeofGenotype(p1,p2)
            g2_alleletype=typeofGenotype(u1,u2)
            #print g1[i], g1_alleletype
            #print g2[i], g2_alleletype
            compared_imputed_genotypes.append( (g1_alleletype, g2_alleletype, g1[i][0])  )
    return  compared_imputed_genotypes




def compare_imputed_genotypes( g1, g2, formatstr, thresh):
    """ given two lists of the form [ (sample, gt_string), ... ] iterate and compare *imputed* genotypes in g1 to genotypes in g2  """
    """ to find out which genotypes in g1 are imputed the FORMAT string must contain the OG tag and the OG must be a nocall """
    """ return list of [ ( g1_alleletype, g2_alleletype, samplename) ] where alleletype is [0, homref; 1, het; 2 hom_nonref; 3, nocall  """
    
    formatfields=formatstr.split(':')
    if 'GPROB' in formatstr and 'OG' in formatstr:
        gprobs_index=formatfields.index('GPROB')
        og_index=formatfields.index('OG')
    else:
        sys.stderr.write("genotype format doesn't contain GPROB/OG, cannot determine if genotype was imputed!")
        exit(1)
    compared_imputed_genotypes=[]

    if len(g1) != len(g2):
        sys.stderr.write("cannot compare phased unphased genotypes; unequal size of lists!")
    else:
        for i in range( 0, len( g1 ) ):

            if g1[i][0] != g2[i][0]:
                sys.stderr.write("samples don't match in genotypes comparison!")
                exit(1)

            #check if the original genotype in g1 was imputed; if not pass on the comparison!
            (og1, og2) = returnAlleles ( getFormatfield( g1[i][1], og_index ) )
            if (typeofGenotype(og1,og2) != 3):
                continue
            g1_maxprob= max ( getFormatfield( g1[i][1], gprobs_index ).split(';') )
            g1_maxprob=float(g1_maxprob)
            #posterior prob of imputed genotype did not meet threshold!
            if g1_maxprob <=float(thresh):
                continue
            (p1,p2) = returnAlleles ( stripGT( g1[i][1] ) )
            (u1, u2) = returnAlleles ( stripGT( g2[i][1]) )

            g1_alleletype=typeofGenotype(p1,p2)
            g2_alleletype=typeofGenotype(u1,u2)
            
            compared_imputed_genotypes.append( (g1_alleletype, g2_alleletype, g1[i][0])  )

    
    return  compared_imputed_genotypes

def compare_genotypes(g1, g2):
    """ given two lists of the form [ (sample, gt_string), ....  ], iterate thru and compare genotypes per individual    """
    """ return list of [ ( g1_alleletype, g2_alleletype, samplename) ] where alleletype is [0, homref; 1, het; 2 hom_nonref; 3, nocall """


    
    compared_genotypes=[] # [ (g1_type, g2_type, sample) ... ]

    if len(g1) != len(g2):
        sys.stderr.write("cannot compare phased unphased genotypes; unequal size of lists!")
    else:
        for i in range( 0, len( g1 ) ):
            if g1[i][0] != g2[i][0]: 
                sys.stderr.write("samples don't match in genotypes comparison!")
                exit(1)
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


def doPatternSearch(pattern, string):
    """ given a pattern object that represents a regular expression search for the pattern in string
        return matched string. if there is no match return None                                  """

    match=pattern.search(string)

    if match.group() == None:
        return None
    else:
        return match.group()