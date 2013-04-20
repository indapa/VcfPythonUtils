#!/usr/bin/env python
from itertools import *
from VcfFile import *
import numpy as np
import re, os
from optparse import OptionParser
from common import grouper, typeofGenotype, melt_lol


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
    #print samples
    #print "#setname\t" + "\t".join(samples) 
    for vrec in vcfobj.yieldVcfRecordwithGenotypes(vcfh):
        outputline=[[ vrec.getPos()] ]
        
        setname=vrec.returnInfoDict()['set'] #which callset does the site belong to?
        
        outputline.append(  [setname] ) #we aggregate genotypes per sample heere
        
        vrec_ziptuple=vrec.zipGenotypes(samples)
        #print vrec_ziptuple
        for (compare, eval) in grouper(2,vrec_ziptuple):
            (comp_allele1, comp_allele2)=compare[1].getAlleles()
            (eval_allele1, eval_allele2)=eval[1].getAlleles()
            eval_alleletype=typeofGenotype(eval_allele1, eval_allele2)
            comp_alleletype=typeofGenotype(comp_allele1, comp_allele2)
            if eval_alleletype == comp_alleletype:
                continue
            outputline.append(   [ eval[0], str(eval_alleletype), compare[0], str(comp_alleletype) ] )
        
        
        print "\t".join (  melt_lol(outputline) )
            
        
        

if __name__ == "__main__":
    main()