#!/usr/bin/python
import sys
import os
import string
import re
from optparse import OptionParser
from VcfFile import *

def main():
    usage = "usage: %prog [options] file.vcf"
    parser = OptionParser(usage)
    (options, args)=parser.parse_args()

    vcfilename=args[0]
    vcfh=open(vcfilename,'r')

    #instantiate a VcfFile object
    vcfobj=VcfFile(vcfilename)
    #parse its metainfo lines (ones that begin with ##)
    vcfobj.parseMetaLines(vcfh)

    #print meta INFO, FORMAT, and FILTER lines
    #vcfobj.addMetaInfoHeader("TV", "Integer", 0, "dummy info")
    vcfobj.printMetaLines()
    

    vcfh.seek(0)

    #parse the header  line #CHROM and print it
    vcfobj.parseHeaderLine(vcfh)
    vcfobj.printHeaderLine()
    samplelist = vcfobj.getSampleList()
    
    
    for vrec in vcfobj.yieldVcfRecordwithGenotypes(vcfh ):
        vrec.checkInfoIds(vcfobj.getMetaInfoNames())
        vrec.check_genotypeFormat(vcfobj.getMetaFormatNames() )
        genotypes = vrec.getGenotypesAlleles()
        sample_genotypelist = zip(samplelist, genotypes)
        callrate=vrec.siteCallrate()
        print str(callrate)
        if vrec.allHets() == False:
            print vrec.toString()
            print sample_genotypelist
        #print sample_genotypelist
        #print vrec.toStringwithGenotypes()
    #print datalines with genotypes included
    vcfobj.printDataLineWithGenotypes(vcfh)

    

if __name__ == "__main__":
    main()
