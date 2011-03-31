#!/usr/bin/python
import sys
import os
import string
import re
from optparse import OptionParser
from vcfIO import *

def main():
    usage = "usage: %prog [options] file.vcf\nConvert VCF to a minimum bed3 format\n"
    parser = OptionParser(usage)
    #parser.add_option("--name", type="string|int|boolean", dest="name", help="help string")

    (options, args)=parser.parse_args()
    vcfile=args[0]

    try:
        vcf_fh=open(vcfile, 'r')
    except:
        sys.stderr.write("unable to open vcf file!\n")

    for vcfline in get_vcfdataline_passfilter(vcf_fh):
        vcf_data=split_vcfdataline(vcfline)
        (chr,pos)=vcf_data[0:2]
        pos=int(pos)
        outstring ="\t".join([ chr, str(pos-1), str(pos) ] )
        print outstring


if __name__ == "__main__":
    main()
