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
    parser.add_option("--addchr", action="store_true", dest="addchr",  help="pre-pend 'chr' to chrom column ", default=False)
    parser.add_option("--dump", action="store_true", dest="dump", help="dump everything after teh ID column in the 4th bed column")
    (options, args)=parser.parse_args()

    vcfilename=args[0]
    #basename, extension = os.path.splitext(vcfilename)
    #bedfile=basename+".bed"
    #bedfh=open(bedfile,'w')
    vcfh=open(vcfilename,'r')
    #instantiate a VcfFile object
    vcfobj=VcfFile(vcfilename)
    #parse its metainfo lines (ones that begin with ##)
    vcfobj.parseMetaAndHeaderLines(vcfh)

    for dataline in vcfobj.yieldVcfDataLine(vcfh):
        fields=dataline.strip().split('\t')
        (chrom,pos,id,ref,alt,qual,filtercode,info)=fields[0:8]
        
        if options.addchr ==True:
            chrom='chr'+chrom
        if filtercode != options.filter and options.filter != None : continue
        (start,end) = (int(pos)-1, int(pos))
        if options.dump == True:
            gstrings=",".join(fields[8::])
            dumpstring=";".join([ref,alt,qual,filtercode,info,gstrings])
            bedstring= "\t".join( [ chrom, str(start), str(end), id ,dumpstring] )
        else:
            bedstring= "\t".join( [ chrom, str(start), str(end), id] )

        print bedstring
if __name__ == "__main__":
    main()
