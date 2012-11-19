#!/usr/bin/env python
import sys
import os
import string
import re
from optparse import OptionParser
from VcfFile import *


def main():
    """ extract homozygous reference records from a VCF file  """
    usage = "usage: %prog [options] file.vcf"
    parser = OptionParser(usage)
    #parser.add_option("--name", type="string|int|boolean", dest="name", help="help string")

    (options, args)=parser.parse_args()
    vcfilename=args[0]
    vcfobj=VcfFile(vcfilename)
    vcfh=open(vcfilename,'r')
    vcfobj.parseMetaAndHeaderLines(vcfh)

    for vrec in vcfobj.yieldVcfRecordwithGenotypes(vcfh):
        if vrec.getAlt() == ".":
            print vrec.toStringwithGenotypes()
        


if __name__ == "__main__":
    main()
