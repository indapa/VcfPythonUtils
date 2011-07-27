#!/usr/bin/env python
import sys
import os
import string
import re
from optparse import OptionParser
from VcfFile import *

""" print the nuumber of each type of variant class ( e.g. snp insertion, deletion, mnp, complex in a VCF file """

def main():
    usage = "usage: %prog [options] arg"
    parser = OptionParser(usage)
    parser.add_option("--info", type="string", dest="infotag", help="INFO tag id that annotates what type of variant the VCF record is", default="TYPE")
    parser.add_option("--filter", type="string", dest="filter", help="only analyze records with matching filter (default is PASS)", default="PASS")

    (options, args)=parser.parse_args()
    if options.infotag == "":
        sys.stderr.write("provide a value for --info parameter!\n")
        exit(1)


    variant_dict={}

    vcfilename=args[0]
    vcfh=open(vcfilename,'r')

    #instantiate a VcfFile object
    vcfobj=VcfFile(vcfilename)
    #parse its metainfo lines (ones that begin with ##)
    vcfobj.parseMetaLines(vcfh)
    descriptors = vcfobj.getMetaInfoDescription()
    infoids=[]
    for (tag, description) in descriptors:
        infoids.append(tag)

    if options.infotag  not in infoids and options.infotag != 'QUAL':
        sys.stderr.write(options.infotag + " tag not in ##INFO headers!\n")
        exit(1)

    vcfh.seek(0)
    vcfobj.parseHeaderLine(vcfh)

    pattern=options.infotag+'=(\w+)'

    for vrec in vcfobj.yieldVcfRecord(vcfh):
        if vrec.getFilter() != options.filter: continue
        searchresult=re.search(pattern, vrec.getInfo() )
        if re.search(pattern, vrec.getInfo() ) == None:
            continue
        else:
            value=re.search(pattern, vrec.getInfo() ).groups()[0]
            if value not in variant_dict.keys():
                variant_dict[value]=1
            else:
                variant_dict[value]+=1

    sys.stderr.write("types and count of different variant classes in " + vcfilename + "\n")
    for k in variant_dict.keys():
        print k, variant_dict[k]
if __name__ == "__main__":
    main()
