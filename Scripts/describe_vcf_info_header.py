#!/usr/bin/python
import sys
import os
import string
import re
from VcfFile import *
from optparse import OptionParser


def main():

    """ prints the description of ##INFO metalines in a VCF  """
    
    usage = "usage: %prog [options] file.vcf"
    parser = OptionParser(usage)
    parser.add_option("--infotag", type="string", dest="infotag", help="prints the  description for the INFO id infotag")
    parser.add_option("--all", action="store_true", dest="all",  help="prints  the  description for  *every* INFO  tag in VCF")
    parser.add_option
    (options, args)=parser.parse_args()
    
    vcfilename=args[0]
    vcfh=open(vcfilename, 'r')
    #instantiate a VcfFile object
    vcfobj=VcfFile(vcfilename)

    #parse its metainfo lines (ones that begin with ##)
    vcfobj.parseMetaLines(vcfh)

    descriptors = vcfobj.getMetaInfoDescription()
    found_tag=0
    for (id, description) in descriptors:
        if options.all==True:
            print id, description
            found_tag=1
            continue
        if id == options.infotag:
            print id, "\t", description
            found_tag=1
    if found_tag  ==0  : sys.stderr.write(options.infotag + " not in ##INFO headers\n")

if __name__ == "__main__":
    main()
