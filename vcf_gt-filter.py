#!/usr/bin/env python

import sys
from optparse import OptionParser
from collections import defaultdict
from VcfFile import *
import argparse

""" filter records to match genotypes according to a pattern
Specify genotypes with the -gt option: -gt <sample><single-space><genotype string>
i.e. -gt "sampleOne 0/0 """

def main():
    usage = "usage: %prog [options] file.vcf "
    parser = argparse.ArgumentParser(description='filter records based on genotypes')
   
    parser.add_argument('vcf', metavar='vcf', type=str,
                   help='vcf file')
    """ http://stackoverflow.com/a/15008806/1735942 """
    parser.add_argument('--no-header',dest='header',action='store_false')
    parser.add_argument('-gt', metavar='gt', type=str, nargs='*', action='append',
                   help='sample 0/0')
   
    args = parser.parse_args()
    
    """ http://stackoverflow.com/q/12460989/1735942 """
    args.gt = [el for elements in args.gt for el in elements]
    
    #print args.gq
    
    gt_filter=[ tuple(x.split(' ')) for x in args.gt ]
    
    gt_dict=defaultdict(list)
    for (k,v) in gt_filter:
        gt_dict[k].append(v)
        
    #print gt_dict
    
    
    
    
    vcfh=open(args.vcf,'r')
    vcfobj=VcfFile(args.vcf)
    vcfobj.parseMetaAndHeaderLines(vcfh)
    header=vcfobj.returnHeader()
    if args.header == True:
        print header
    samplelist=vcfobj.getSampleList()   
    for s in gt_dict.keys():
        if s not in samplelist:
            print s ," not in samples!\n"
            sys.exit(1)
    #print header
    #print header
    #print gt_dict.keys()

    for vrec in vcfobj.yieldVcfRecordwithGenotypes(vcfh ):
        genotypes_toFilter=[] #list [ (sample,genoObj) ] to be filtered 
        genotype_tuple= vrec.zipGenotypes(samplelist) ## get a list of tuples [ (sample, VcfGenotype object) ... ]
        for (s,g) in genotype_tuple:
            if s in  gt_dict.keys():
                #print s
                if  len(gt_dict[s]) > 1: # logical or 
                    
                    if any( [ g.getFormatVal('GT') == v for v in gt_dict[s] ] ):
                        genotypes_toFilter.append(True)
                    else: genotypes_toFilter.append(False)
                else:
                    if all( [ g.getFormatVal('GT') == v for v in gt_dict[s] ] ):
                        genotypes_toFilter.append(True)
                    else:genotypes_toFilter.append(False)
                
        #print genotypes_toFilter
        if all(item == True for item in genotypes_toFilter):
            print vrec.toStringwithGenotypes()
                
        
        

if __name__ == "__main__":
    main()
