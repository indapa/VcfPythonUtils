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
    parser.add_option("--v", action="store_true", dest="snp",  help="restrict analysis to SNPs (must have INFO ID SNP in header")

    (options, args)=parser.parse_args()

    vcfilename=args[0]
    vcfh=open(vcfilename,'r')

    #instantiate a VcfFile object
    vcfobj=VcfFile(vcfilename)
    #parse its metainfo lines (ones that begin with ##)
    vcfobj.parseMetaLines(vcfh)
    vcfh.seek(0)
    vcfobj.parseHeaderLine(vcfh)

    for vrec in vcfobj.yieldVcfRecord(vcfh)
        print vrec.getInfo()


if __name__ == "__main__":
    main()
