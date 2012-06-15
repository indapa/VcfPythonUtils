#!/usr/bin/env python
import sys
import os
import string
import re
from optparse import OptionParser
from VcfFile import *
from collections import *

""" given a vcf that has been merged with GATK CombineVariants, get the counts of variants in each portion of the venn; with the --extract option print the records belonging to a particualr set """
def main():
    usage = "usage: %prog [options] merged.vcf"
    parser = OptionParser(usage)
    parser.add_option("--extract", type="string", dest="extract", help="set to extract", default=None)

    (options, args)=parser.parse_args()
    if options.extract != None:
        pattern='set=('+options.extract+')'
    else:
        pattern='set=(\w+)'
    
    callsetCounter=Counter()

    vcfobj=VcfFile('gatkuf.20120527.sample3_4.wga.genomic.merge.vcf')
    vcfh=open('gatkuf.20120527.sample3_4.wga.genomic.merge.vcf','r')

    vcfobj.parseMetaAndHeaderLines(vcfh)
    if options.extract != None:
        vcfobj.printMetaAndHeaderLines()
    samples=vcfobj.getSampleList()
    for vrec in vcfobj.yieldVcfRecordwithGenotypes(vcfh):
    
        if re.search(pattern, vrec.getInfo() ) == None:
            continue
        else:
            value=re.search(pattern, vrec.getInfo() ).groups()[0]
        #print value
            callsetCounter[value]+=1
        if options.extract != None:
            print vrec.toStringwithGenotypes()

    if options.extract == None:
        for (type,count) in callsetCounter.items():
            print type, str(count)
        print "total: " + str( sum(callsetCounter.values()) )



if __name__ == "__main__":
    main()
