#!/usr/bin/python
import sys
import os
import string
import re
from optparse import OptionParser

from VcfFile import *

""" look for shared variants in a list of affecteds given an inheritance model """

def main():
    usage = "usage: %prog [options] file.vcf"
    parser = OptionParser(usage)
    parser.add_option("--model", type="string", dest="model", default = "dominant", help=" inheritance model (dominant recessive) ")
    parser.add_option("--affected", type="string", dest="affected", help="sample name of affecteds (one per line)")
    parser.add_option("--filter", type="string", dest="filter", help="analyze only those  records matching filter (default is None)", default=None)
    (options, args)=parser.parse_args()
    if options.affected==None:
        sys.stderr.write("please provide a value to --affected parameter!\n")
        exit(1)

    affecteds=[] # list of affected samples
    affectedfh=open(options.affected, 'r')
    for line in affectedfh:
        affecteds.append(line.strip() )

    vcfilename=args[0]
    vcfh=open(vcfilename,'r')

    #instantiate a VcfFile object
    vcfobj=VcfFile(vcfilename)
    samplelist=vcfobj.getSampleList()
    #parse its metainfo lines (ones that begin with ##)
    vcfobj.parseMetaLines(vcfh)

    for vrec in vcfobj.yieldVcfRecordwithGenotypes(vcfh ):
        if vrec.getFilter() != options.filter and options.filter != None : continue
        genotype_tuple= vrec.zipGenotypes(samplelist)
        print genotype_tuple

if __name__ == "__main__":
    main()
