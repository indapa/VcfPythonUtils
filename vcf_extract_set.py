#!/usr/bin/env python
import gzip
import sys,os,re
from common import * 

from VcfFile import *
import argparse
import os


def main():
    
    """  This program extracts out records matching the set=(\S+). Typically, the VCF is derived
    from GATK CombineVariants, but any vcf with set=(\S+) can be examined with this program """
   
    usage = "usage: %prog [options] file.vcf.gz "
    #parser = OptionParser(usage)
    parser = argparse.ArgumentParser(description=' extract records with matching set=(\S+) tag')
    
    parser.add_argument('vcfile', metavar='vcfile', type=str, help='file.vcf.gz')
    #parser.add_argument('-filter', dest='filter', type=str, default=".", help='filter value')
    parser.add_argument('-set', dest='set', type=str, default=None, help="name of set to extract")
    
    args = parser.parse_args()
    if args.set == None: 
        sys.stderr.write("please provide value to -set option!\n")
        sys.exit(1)
   
    (path, vcfile)=os.path.split(args.vcfile )
    
    basename=return_file_basename( return_file_basename(vcfile) )
    sys.stderr.write( basename +"\n")
    
    outvcf=".".join([basename, args.set, 'vcf'])
    sys.stderr.write( outvcf +"\n")
    outfh=open(outvcf, 'w')
    
    if args.vcfile.endswith(".gz"):
        vcfh=gzip.open(args.vcfile,'r')
    else:
        vcfh=open(args.vcfile,'r')
    vcfobj=VcfFile(args.vcfile)
    
    pattern=';set=(\S+)'
   
    vcfobj.parseMetaAndHeaderLines(vcfh)
    header=vcfobj.returnHeader() 
    outfh.write( header +"\n")
    
    for vrec in vcfobj.yieldVcfRecordwithGenotypes(vcfh):
        searchresult=re.search(pattern, vrec.getInfo() )
        if re.search(pattern, vrec.getInfo() ) == None:
            continue
        value=re.search(pattern, vrec.getInfo() ).groups()[0]
        #print value
        if value == args.set:
            outfh.write(  vrec.toStringwithGenotypes() +"\n" )    
    

if __name__ == "__main__":
    main()
