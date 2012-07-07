#!/usr/bin/env python
import sys
import os
import string
import re
from optparse import OptionParser
import collections
from VcfFile import *

import numpy as np
from common import *



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

    TsTv_counter=collections.Counter()
    RefAlt_counter=collections.Counter()
   

    



    samples=vcfobj.getSampleList()


    genotype_dict={}
    for s in samples:
        genotype_dict[s]=[0,0,0,0]

    for vrec in vcfobj.yieldVcfRecordwithGenotypes(vcfh):
        if vrec.getFilter() != options.filter and options.filter != None:
            
            continue
        #print vrec.toString()


        ref=vrec.getRef()
        numAlleles=vrec.getAlt().split(',')
       
        for alt in numAlleles:
            if isTransition(ref,alt) == True:
                TsTv_counter['transition']+=1
            else:
                TsTv_counter['transversion']+=1
            refalt_string=" ".join( [ ref, alt])
            if len(alt) ==1 and len(ref) ==1:
                RefAlt_counter[ refalt_string ]+=1
        #    sys.stderr.write("need to add code to accomodate non-biallelic sites... skipping record\n")
        #    continue
        vrec_ziptuple=vrec.zipGenotypes(samples)
        genotype_typecounts=get_genotype_counts(vrec_ziptuple)
        for (g, sample) in genotype_typecounts:
            #print g,sample
            if g == None: continue
            genotype_dict[sample][g]+=1
    print 
    print "#sample homoz_ref het homoz_nonref nocall"
    for sample in genotype_dict.keys():
        
        outstring = " ".join( map(str,genotype_dict[sample]) )
        print sample, outstring 

    print

    for (type,count) in TsTv_counter.items():
        print type, count
    print sum(TsTv_counter.values())
    print
    for (type, count) in RefAlt_counter.items():
        print type, count
    print sum(RefAlt_counter.values())
if __name__ == "__main__":
    main()
