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
from itertools import *


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
    counter=0
    for vrec in vcfobj.yieldVcfRecordwithGenotypes(vcfh):
        if vrec.getFilter() != options.filter and options.filter != None:
            sys.stderr.write("skipped filter..\n")
            continue
        #print vrec.toString()
        counter+=1
        
        if vrec.getAlt() == ".": continue
        ref=vrec.getRef()
        numAlleles=vrec.getAlt().split(',')
        
        if len(numAlleles) > 1:
            sys.stderr.write("multi alleleic record\n")
         
       
        for alt in numAlleles:
            if len(alt) ==1 and len(ref) ==1:
                if isTransition(ref,alt) == True:
                    TsTv_counter['transition']+=1
                else:
                    TsTv_counter['transversion']+=1
                refalt_string=" ".join( [ ref, alt])
            #since the number of alleles on indels is unbounded, we only keep track of single nucleotide substitutions
                RefAlt_counter[ refalt_string ]+=1
                
        
        vrec_ziptuple=vrec.zipGenotypes(samples)
        genotype_typecounts=get_genotype_counts(vrec_ziptuple)
        for (g, sample) in genotype_typecounts:
            #print g,sample
            if g == None:
                sys.stderr.write("skipped genotype\n")
                continue
            genotype_dict[sample][g]+=1

    
    print
    print " ".join( ['sample', 'homoz_ref', 'het', 'homoz_nonref', 'nocall', 'total'])
    for sample in genotype_dict.keys():
        """ http://docs.python.org/library/functions.html#reduce """
        tota=reduce(lambda x, y: x+y,genotype_dict[sample])
        
        outstring = " ".join( map(str,genotype_dict[sample]) )
        print " ".join ( [sample, outstring,str(tota)])

    print

    for (type,count) in TsTv_counter.items():
        print type, count
    
    

    totalpercent=0
    for a1,a2  in combinations('ACGT',2):
        count1 = RefAlt_counter[ ' '.join ( [ a1, a2] ) ]
        count2 = RefAlt_counter[ ' '.join ( [ a2, a1] ) ]
        total=count1 + count2
        try:
            percent= round ( float(total) / float(sum(RefAlt_counter.values()) ), 4)
            print ' '.join ( [ a1, a2] ), str(total), str(percent)
            totalpercent+=percent
        except ZeroDivisionError:
            sys.stderr.write( " integer division or modulo by zero\n")
    #for (type, count) in RefAlt_counter.items():
    #    print type, count
    print sum(RefAlt_counter.values()), str(totalpercent)


    print "Total vcf records: " + str(counter) + "\n"
if __name__ == "__main__":
    main()
