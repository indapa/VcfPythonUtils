#!/usr/bin/python
import sys
import os
import string
import re
from optparse import OptionParser
from VcfFile import *
from collections import *
from itertools import *
import numpy as np


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
    usage = "usage: %prog [options] nrd.log.vcf"
    parser = OptionParser(usage)
    (options, args)=parser.parse_args()
    vcfile= args[0]
    vcfobj=VcfFile(vcfile)
    vcfh=open(vcfile, 'r')

    vcfobj.parseMetaAndHeaderLines(vcfh)
    header=vcfobj.returnHeader()
    #print header
    samples=vcfobj.getSampleList()


    altallele_counter=Counter()

    nrdtable= np.matrix( [ [ 0,0 ], [ 0,0 ] ] )
    print "\t".join( [ 'chrom', 'pos', 'evalref', 'evalalt', 'evalalt_perct', 'compref', 'compalt', 'compalt_perct','evalgt', 'compgt'])
    for vrec in vcfobj.yieldVcfRecordwithGenotypes(vcfh):
        if len(vrec.getAlt()) > 1: continue #need to fix this to deal with non-biallelic sites
        if 'set=Intersection' not in vrec.getInfo(): continue # all sites should be in intersection
        vrec_ziptuple=vrec.zipGenotypes(samples)

        for (compare, eval) in grouper(2,vrec_ziptuple):
            (comp_allele1, comp_allele2)=compare[1].getAlleles()
            (eval_allele1, eval_allele2)=eval[1].getAlleles()
            eval_alleletype=typeofGenotype(eval_allele1, eval_allele2)-1
            comp_alleletype=typeofGenotype(comp_allele1, comp_allele2)-1


            comp_ad=compare[1].getFormatVal('AD')
            eval_ad=eval[1].getFormatVal('AD')

            

            nrdtable[eval_alleletype, comp_alleletype]+=1

            if eval_alleletype == 0:
                if comp_alleletype == 1:
                    #print vrec.toStringwithGenotypes()
                    #print comp_ad, eval_ad
                    altallele_counter[ "\t".join( [ vrec.getRef(), vrec.getAlt() ] ) ]  +=1
                    (comp_ref, comp_alt)=comp_ad.split(',')
                    (eval_ref, eval_alt)=eval_ad.split(',')

                    comp_altpercentage= float(comp_alt)/ (float(comp_alt)+float(comp_ref))
                    comp_altpercentage=round(comp_altpercentage,2)

                    eval_altpercentage= float(eval_alt)/ (float(eval_alt)+float(eval_ref))
                    eval_altpercentage=round(eval_altpercentage,2)



                    outstring="\t".join([ vrec.getChrom(), vrec.getPos(), eval_ref, eval_alt, str(eval_altpercentage),  comp_ref, comp_alt, str(comp_altpercentage), str( typeofGenotype(eval_allele1, eval_allele2) ), str( typeofGenotype(comp_allele1, comp_allele2) ) ])
                    print outstring
    #print nrdtable

    for (type,count) in altallele_counter.items():
        print type,count

if __name__ == "__main__":
    main()
