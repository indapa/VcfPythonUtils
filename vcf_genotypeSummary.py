#!/usr/bin/env python
import sys
import os
import string
import re
from optparse import OptionParser

from VcfFile import *
import numpy as np

def typeofGenotype(allele1, allele2):
    """ I really should be a python version of a typedef here, but dont know how
        hom_ref =0 het =1 hom_nonref=2 no_call=3 """

    if allele1 == '0' and allele2 == '0':
        return 0
    elif allele1 == '0' and allele2== '1':
        return 1
    elif allele1 =='1' and allele2 == '0':
        return 1
    elif allele1== '1' and allele2== '1':
        return 2
    elif allele1== '.' or allele2 == '.':
        return 3
    elif allele1 == '.' and allele2 == '.':
        return 3
    else:
        return None


def get_genotype_counts(g):
    genotype_counts=[]
    for i in range(0, len(g)):
        (p1,p2)=g[i][1].getAlleles()
        alleletype=typeofGenotype(p1,p2)
        genotype_counts.append( (alleletype, g[i][0] ) )

    return genotype_counts

def main():
    usage = "usage: %prog [options] file.vcf"
    parser = OptionParser(usage)
    parser.add_option("--filter", type="string", dest="filter", help="analyze only those records with matching filter")

    (options, args)=parser.parse_args()

    vcfilename=args[0]
    vcfh=open(vcfilename,'r')

    #instantiate a VcfFile object
    vcfobj=VcfFile(vcfilename)
    #parse its metainfo lines (ones that begin with ##)
    vcfobj.parseMetaAndHeaderLines(vcfh)

    samples=vcfobj.getSampleList()


    genotype_dict={}
    for s in samples:
        genotype_dict[s]=[0,0,0,0]

    for vrec in vcfobj.yieldVcfRecordwithGenotypes(vcfh):
        if vrec.getFilter() != options.filter and options.filter != None:
            continue
        #print vrec.toString()
        numAlleles=vrec.getAlt().split(',')
        if len(numAlleles) > 1:
            sys.stderr.write("need to add code to accomodate non-biallelic sites... skipping record\n")
            continue
        vrec_ziptuple=vrec.zipGenotypes(samples)
        genotype_typecounts=get_genotype_counts(vrec_ziptuple)
        for (g, sample) in genotype_typecounts:
            #print g,sample
            if g == None: continue
            genotype_dict[sample][g]+=1
    print "#sample homoz_ref het homoz_nonref nocall"
    for sample in genotype_dict.keys():
        
        outstring = " ".join( map(str,genotype_dict[sample]) )
        print sample, outstring 


if __name__ == "__main__":
    main()
