#!/usr/bin/env python
from itertools import *
from VcfFile import *
import numpy as np
import re, os
from optparse import OptionParser
from common import grouper
from common import typeofGenotype

def main():
    usage = "usage: %prog [options]  nrd.log.vcf\n" 
    parser = OptionParser(usage)
    #parser.add_option("--matrixonly", action="store_true", dest="matrixonly", help="only print concordance matrixe", default=False)
    #parser.add_option("--includeRef", action="store_true", dest="includeRef", help="include sites in the set ReferenceInAll", default=False)

    (options, args)=parser.parse_args()
    vcfilename=args[0]
    basename=os.path.splitext(vcfilename)[0]
    
    vcfobj=VcfFile(vcfilename)
    vcfh=open(vcfilename,'r')

    vcfobj.parseMetaAndHeaderLines(vcfh)
    samples=vcfobj.getSampleList()
    print samples
    
    for vrec in vcfobj.yieldVcfRecordwithGenotypes(vcfh):
        vrec_ziptuple=vrec.zipGenotypes(samples)
        for (compare, eval) in grouper(2,vrec_ziptuple):
            print compare
            print eval
            print
        print
        

if __name__ == "__main__":
    main()