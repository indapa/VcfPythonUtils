#!/usr/bin/env python
import sys
import os
import string
import re
from optparse import OptionParser

from VcfFile import *

""" print records belonging to a certain type of variant class (e.g. SNP) in a VCF file """

def main():
    usage = "usage: %prog [options] file.vcf\n print records belonging to a certain type of variant class (e.g. SNP) in a VCF file\n\n"
    parser = OptionParser(usage)
    parser.add_option("--info", type="string", dest="infotag", help="INFO tag id that annotates what type of variant the VCF record is", default="TYPE")
    parser.add_option("--type", type="string", dest="variantype", help="type of variant (SNP INS DEL)", default="")
    parser.add_option("--filter", type="string", dest="filter", help="extract records matching filter (default is None)", default=None)

    (options, args)=parser.parse_args()
    if options.infotag == "":
        sys.stderr.write("provide a value for --info parameter!\n")
        exit(1)
    if options.variantype == "":
        sys.stderr.write("provide a value of --type parameter!\n")
        exit(1)

    variant_dict={}

    vcfilename=args[0]
    vcfh=open(vcfilename,'r')

    #instantiate a VcfFile object
    vcfobj=VcfFile(vcfilename)
    #parse its metainfo lines (ones that begin with ##)
    vcfobj.parseMetaLines(vcfh)
    vcfobj.printMetaLines()
    descriptors = vcfobj.getMetaInfoDescription()
    infoids=[]
    for (tag, description) in descriptors:
        infoids.append(tag)

    if options.infotag  not in infoids and options.infotag != 'QUAL':
        sys.stderr.write(options.infotag + " tag not in ##INFO headers!\n")
        exit(1)

    vcfh.seek(0)
    vcfobj.parseHeaderLine(vcfh)
    vcfobj.printHeaderLine()

    pattern=options.infotag+'=('+options.variantype+')'

    for vrec in vcfobj.yieldVcfRecord(vcfh):
        if vrec.getFilter() != options.filter and options.filter != None : continue
        searchresult=re.search(pattern, vrec.getInfo() )
        if re.search(pattern, vrec.getInfo() ) == None:
            continue
        else:
            value=re.search(pattern, vrec.getInfo() ).groups()[0]
            print vrec.toString()
if __name__ == "__main__":
    main()
