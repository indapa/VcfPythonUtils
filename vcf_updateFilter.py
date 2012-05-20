#!/usr/bin/python
import sys
import os
import string
import re
from optparse import OptionParser

from VcfFile import *

""" update all records in vcf file with filter value specified by --filter to the one specified by --newfilter  """
def main():
    usage = "usage: %prog [options] file.vcf\n print records belonging to a certain type of variant class (e.g. SNP) in a VCF file\n\n"
    parser = OptionParser(usage)
    parser.add_option("--info", type="string", dest="infotag", help="INFO tag id that annotates what type of variant the VCF record is", default="TYPE")

    parser.add_option("--filter", type="string", dest="filter", help="extract records matching filter (default is None)", default=None)
    parser.add_option("--newfilter", type="string", dest="newfilter", help="new value for filter column", default=None)
    
    (options, args)=parser.parse_args()
    if options.filter == None:
        sys.stderr.write("provide a value of --filter parameter!\n")
        exit(1)
    if options.newfilter == None:
        sys.stderr.write("provide a value of --newfilter parameter!\n")
        exit(1)

    

    vcfilename=args[0]
    vcfh=open(vcfilename,'r')

    #instantiate a VcfFile object
    vcfobj=VcfFile(vcfilename)
    #parse its metainfo lines (ones that begin with ##)
    vcfobj.parseMetaAndHeaderLines(vcfh)
    vcfobj.printMetaAndHeaderLines()

    for vrec in vcfobj.yieldVcfRecordwithGenotypes(vcfh):
        if vrec.getFilter() == options.filter:
            vrec.setFilter(options.newfilter)
            
        print vrec.toStringwithGenotypes()


if __name__ == "__main__":
    main()
