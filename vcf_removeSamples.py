#!/usr/bin/env python
import gzip
from itertools import *
from VcfFile import *
from VcfSampleEval import *
from optparse import OptionParser
import argparse
import os


def main():
    
    """  remove samples from a vcf file """
    usage = "usage: %prog [options] file.vcf.gz "
    #parser = OptionParser(usage)
    parser = argparse.ArgumentParser(description='remove samples from vcf file')
    parser.add_argument('removesamples', metavar='sample', type=str, nargs='+',
                   help='sample names to remove')
    parser.add_argument('-vcf', dest='vcfile', type=str, help="vcf file to remove samples from")
    #parser.add_argument("vcf", help="vcf file to analyze")
    args = parser.parse_args()
    #print 'remove these samples: ', args.samples
    #print args.vcfile
    
    
    vcfh=gzip.open(args.vcfile,'r')
    vcfobj=VcfFile(args.vcfile)
    
    vcfobj.parseMetaAndHeaderLines(vcfh)
    
    #print header
    samples=vcfobj.getSampleList()
    newsamples= [ s for s in samples if s not in args.removesamples]
    #print 'keep these samples: ',  newsamples
    vcfobj.setSampleList(newsamples)
    header=vcfobj.returnHeader() 
    print header
    
    for vrec in vcfobj.yieldVcfRecordwithGenotypes(vcfh):
        keepGenotypes=[]
        vrec_ziptuple=vrec.zipGenotypes(samples)
        for (s, genObj) in vrec_ziptuple:
            if s not in args.removesamples:
                #print s
                keepGenotypes.append( genObj )
        #print keepGenotypes
        vrec.addGenotypeList(  keepGenotypes )
        print vrec.toStringwithGenotypes()
                
             
        
                    
   



if __name__ == "__main__":
    main()
