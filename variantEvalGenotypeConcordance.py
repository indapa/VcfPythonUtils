#!/usr/bin/env python
from itertools import *
from VcfFile import *
import numpy as np

from optparse import OptionParser

def typeofGenotype(allele1, allele2):
    """ I really should be a python version of a typedef here, but dont know how
        hom_ref =0 het =1 hom_nonref=2 no_call=3                              """

    if allele1 == '0' and allele2 == '0': return 0

    if allele1 == '0' and allele2== '1': return 1
    if allele1 =='1' and allele2 == '0': return 1

    if allele1== '1' and allele2== '1': return 2

    if allele1== '.' or allele2 == '.': return 3
    if allele1 == '.' and allele2 == '.': return 3

""" iterate through an utterable n values at a time
     http://stackoverflow.com/a/2990151         """
def grouper(n, iterable, fillvalue=None):
    "grouper(3, 'ABCDEFG', 'x') --> ABC DEF Gxx"
    args = [iter(iterable)] * n
    return izip_longest(fillvalue=fillvalue, *args)


def main():
    usage = "usage: %prog [options] file.vcf"
    parser = OptionParser(usage)
    (options, args)=parser.parse_args()
    vcfilename=args[0]
#row is veal column is comparison
    concordancetable= np.matrix( [ [ 0,0,0,0 ], [ 0,0,0,0 ], [ 0,0,0,0 ], [ 0,0,0,0 ] ] )



    vcfobj=VcfFile(vcfilename)
    vcfh=open(vcfilename,'r')

    vcfobj.parseMetaAndHeaderLines(vcfh)
    samples=vcfobj.getSampleList()
    for vrec in vcfobj.yieldVcfRecordwithGenotypes(vcfh):
        if len(vrec.getAlt()) > 1: continue
        if 'filterIn' in vrec.getInfo(): 
        #print vrec.toString()
            continue
        vrec_ziptuple=vrec.zipGenotypes(samples)
        for (compare, eval) in grouper(2,vrec_ziptuple):
            (comp_allele1, comp_allele2)=compare[1].getAlleles()
            (eval_allele1, eval_allele2)=eval[1].getAlleles()
            eval_alleletype=typeofGenotype(eval_allele1, eval_allele2)
            comp_alleletype=typeofGenotype(comp_allele1, comp_allele2)
            concordancetable[eval_alleletype, comp_alleletype]+=1
       

    discordance=concordancetable[0,1]+concordancetable[0,2]+concordancetable[1,0]+concordancetable[1,2]+concordancetable[2,0]+concordancetable[2,1]
    total=concordancetable[0,1]+concordancetable[0,2]+concordancetable[1,0]+concordancetable[1,1]+ concordancetable[1,2]+concordancetable[2,0]+concordancetable[2,1] +concordancetable[2,2]

    nrd=round( (float(discordance)/float(total)) * 100, 2)
    variant_count_evaluation= concordancetable[1,1]+ concordancetable[1,2]+ concordancetable[2,1]+ concordancetable[2,2]
    variant_count_comparison= concordancetable[0,1]+concordancetable[0,2]+concordancetable[1,1]+concordancetable[1,2]+concordancetable[2,1]+concordancetable[2,2]+concordancetable[3,1]+concordancetable[3,2]
    nrs=round( float(variant_count_evaluation)/float(variant_count_comparison) * 100 , 2)
    
    print "NRD: ", str(nrd)
    print "NRS ", str(nrs)
# <codecell>

if __name__ == "__main__":
    main()

