#!/usr/bin/python
import sys
import os
import string
import re
from optparse import OptionParser
from VcfFile import *

""" convert a vcf file to a bed4 file *remember that bed files are zero-based, half open """
def main():
    usage = "usage: %prog [options] arg"
    parser = OptionParser(usage)
    parser.add_option("--filter", type="string", dest="filter", help="extract records matching filter (default is None)", default=None)
    parser.add_option("--addchr", action="store_true", dest="addchr",  help="pre-pend 'chr' to output ")
    (options, args)=parser.parse_args()

    vcfilename=args[0]
    basename, extension = os.path.splitext(vcfilename)
    bedfile=basename+".bed"
    bedfh=open(bedfile,'w')
    vcfh=open(vcfilename,'r')
    #instantiate a VcfFile object
    vcfobj=VcfFile(vcfilename)
    #parse its metainfo lines (ones that begin with ##)
    vcfobj.parseMetaLines(vcfh)

    for dataline in vcfobj.yieldVcfDataLine(vcfh):
        fields=dataline.strip().split('\t')
        (chrom,pos,id,ref,alt,qual,filtercode,info)=fields[0:8]
        if options.addchr ==True:
            chrom='chr'+chrom
        if filtercode != options.filter and options.filter != None : continue
        (start,end) = (int(pos)-1, int(pos))
        bedstring= "\t".join( [ chrom, str(start), str(end), id ] )
        bedfh.write(bedstring+"\n")
if __name__ == "__main__":
    main()
