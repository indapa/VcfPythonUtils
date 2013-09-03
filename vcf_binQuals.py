#!/usr/bin/env python
import gzip
from itertools import *
import sys,os
from common import * 

from VcfFile import *
import argparse
import os
import numpy as np





def main():
    
    """  This program bins the records of a VCF file according to a user defined range and number of bins.
        For example if -start 10 and -end 100 and -num of 10 it would make 10 bins:
        10, 20, 30, 40,50,60,70,80,90,100
        
        Then for each record if the QUAL is >=x, then that record is written to *.qual_x.vcf file
    """
    usage = "usage: %prog [options] file.vcf.gz "
    #parser = OptionParser(usage)
    parser = argparse.ArgumentParser(description=' bin vcf records according to QUAL')
    
    parser.add_argument('vcfile', metavar='vcfile', type=str, help='file.vcf.gz')
    parser.add_argument('-filter', dest='filter', type=str, default=".", help='filter value')
    parser.add_argument('-start', dest='start', type=int, help="starting point for QUAL range")
    parser.add_argument('-stop', dest='end', type=int, help="ending pint for QUAL range")
    parser.add_argument('-num', dest='num', type=int, help='number of bins')
    
    #parser.add_argument("vcf", help="vcf file to analyze")
    args = parser.parse_args()
    #print args
    (path, vcfile)=os.path.split(args.vcfile )
    basename=return_file_basename( return_file_basename(vcfile) )
    print basename
    if args.start == None or  args.end  == None or  args.num == None:
        sys.stderr.write("please give start stop and number of bins for QUAL")
        sys.exit(1)
        
    bins=np.linspace(args.start, args.end, args.num)
    binstring=bins
    binstring=binstring.astype(int).tolist()
    print binstring
    binned_vcfilenames=[ ".".join( [ basename, "qual_"+ str(s), "vcf"]) for s in binstring ]
    print binned_vcfilenames
    #binned_fh = itertools.chain(*(open(f, "w") for f in binned_vcfilenames))   
    """ we create a list of filehandles for the binned VCFs """                     
    binned_fh=list(itertools.imap(lambda x:open(x,'w'), binned_vcfilenames))                               
    
    if args.vcfile.endswith(".gz"):
        vcfh=gzip.open(args.vcfile,'r')
    else:
        vcfh=open(args.vcfile,'r')
    vcfobj=VcfFile(args.vcfile)
    
    
    #vcf_reader = vcf.Reader(open(args.vcfile, 'r'))
    #print vcf_reader.metadata
    vcfobj.parseMetaAndHeaderLines(vcfh)
    header=vcfobj.returnHeader() 
    
    map(lambda x: x.write(header+"\n"),binned_fh)
    
    
    #vcfrecord_bins= [ [] for  i in xrange(len(bins)) ]
    sys.stderr.write("binning vcf records based on quality ....\n")
    #for vrec in vcf_reader:
    for vrec in vcfobj.yieldVcfRecordwithGenotypes(vcfh):
        """ skip records that do not have PASS in filter column """
        if vrec.getFilter() != args.filter:
            continue
        QUAL=float(vrec.getQual())
        vcfstring=vrec.toStringwithGenotypes()
        for i in xrange(len(bins)):
            if QUAL >= bins[i]:
                binned_fh[i].write(vcfstring+"\n")
                #vcfrecord_bins[i].append(vrec)
            else: break
        

                
    map(lambda x: x.close(),binned_fh)
    
if __name__ == "__main__":
    main()