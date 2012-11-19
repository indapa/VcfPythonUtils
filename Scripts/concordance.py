#!/usr/bin/env python
import sys
import os
import string
import re
from optparse import OptionParser
from VcfFile import *
def main():
    usage = "usage: %prog [options] arg"
    parser = OptionParser(usage)

    (options, args)=parser.parse_args()

    vcfobj=VcfFile('gatkuf.20120527.sample3_4.wga.genomic.merge.vcf')
    vcfh=open('gatkuf.20120527.sample3_4.wga.genomic.merge.vcf','r')

    vcfobj.parseMetaAndHeaderLines(vcfh)
    for vrec in vcfobj.yieldVcfRecordwithGenotypes(vcfh):
        print vrec.toStringwithGenotypes()


if __name__ == "__main__":
    main()
