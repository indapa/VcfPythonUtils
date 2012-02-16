#!/usr/bin/python
import sys
import os
import string
import re
from optparse import OptionParser

from VcfFile import *

""" remove duplicate records in a vcf file based on chrom and pos of vcf datalines and write duplicate datalines to duplicate.log file"""

def main():
    usage = "usage: %prog [options] file.vcf"
    parser = OptionParser(usage)
    parser.add_option("--filter", type="string", dest="filter", help="only analyze records with matching filter (default is None)", default=None)
    (options, args)=parser.parse_args()


    variant_dict={} #key variant type value VcfRecord object
    vcfilename=args[0]
    vcfh=open(vcfilename,'r')
    dupfh=open('duplicate.log', 'w')
    #instantiate a VcfFile object
    vcfobj=VcfFile(vcfilename)
    #parse its metainfo lines (ones that begin with ##)
    vcfobj.parseMetaLines(vcfh)
    vcfh.seek(0)
    vcfobj.parseHeaderLine(vcfh)

    for dataline in vcfobj.yieldVcfDataLine(vcfh):
        fields=dataline.strip().split('\t')
        (chrom,pos,id,ref,alt,qual,filtercode,info)=fields[0:8]
        if filtercode != options.filter and options.filter != None: continue
      
        stringkey="\t".join([chrom, pos])
        
        if stringkey not in variant_dict.keys():
            variant_dict[stringkey]=1
            print dataline
        else:
           dupfh.write(dataline+"\n")

if __name__ == "__main__":
    main()
