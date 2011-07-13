#!/usr/bin/python
import sys
import os
import string
import re
from optparse import OptionParser

from VcfFile import *

def main():
    usage = "usage: %prog --in file.vcf --in file2.vcf [...] --out out.vcf\n merge multiple vcf files and print to STDOUT"
    vcf_files = []
    parser = OptionParser(usage)
    
    parser.add_option("--in", action="append",  dest="vcf_files", help="outfile")
    (options, args)=parser.parse_args()
    

    commandline = " ".join(sys.argv)
    commandlineheader="##commandline " + commandline 
    print commandlineheader
    if len(options.vcf_files) <= 1:
        sys.stderr.write("need at least two vcf files to merge!\n")
        exit(1)

    printheader=False

    for file in options.vcf_files:
        
        #instantiate a VcfFile object
        vcfh=open(file,'r')
        vcfobj=VcfFile(file)
        #parse its metainfo lines (ones that begin with ##)
        vcfobj.parseMetaLines(vcfh)
        vcfh.seek(0)
        vcfobj.parseHeaderLine(vcfh)

        if printheader ==False:
            vcfobj.printMetaLines()
            vcfobj.printHeaderLine()
            printheader=True
        vcfh.seek(0)

        for dataline in vcfobj.yieldVcfDataLine(vcfh):
            print dataline

        vcfh.close()

if __name__ == "__main__":
    main()
