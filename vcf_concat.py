#!/usr/bin/python
import sys
import os
import string
import re
from optparse import OptionParser

from VcfFile import *

""" concatenate an arbritary number of VCF files into a single file and print to STDOUT
    note - doesn't check for equivalent metadata information across files             """

def main():
    usage = "usage: %prog  --in file.vcf --in file2.vcf [ ... ] "
    parser = OptionParser(usage)
    parser.add_option("--in", action="append", dest="vcf_list", help="file to concatenate")

    (options, args)=parser.parse_args()
    commandline=" ".join(sys.argv)
    print "##commandLine " + commandline
    printmetalines=False
    
    for file in options.vcf_list:
        vcfh=open(file, 'r')
        vcfobj=VcfFile(file)
        #parse its metainfo lines (ones that begin with ##)
        if printmetalines == False:
            vcfobj.parseMetaLines(vcfh)
            vcfobj.printMetaLines()
            vcfh.seek(0)
            vcfobj.parseHeaderLine(vcfh)
            vcfobj.printHeaderLine()
            printmetalines=True
            vcfh.seek(0)
            
        for line in vcfobj.yieldVcfDataLine(vcfh):
            print line

if __name__ == "__main__":
    main()
