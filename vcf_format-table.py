#!/usr/bin/env python
from itertools import *
from VcfFile import *
from VcfSampleEval import *
from optparse import OptionParser
from common import grouper
from common import typeofGenotype
import os

def main():
    usage = "usage: %prog [options] file.vcf \n output format values from  genotype data field  in a VCF  for suitabale plotting/dataviz"
    parser = OptionParser(usage)
    parser.add_option("--includeRef", action="store_true", dest="includeRef", help="include sites in the set ReferenceInAll", default=False)
    parser.add_option("--includeFilter", action="store_true", dest="includeFilter", help="include site filtered or not!", default=False)
    parser.add_option("--formatTag", dest="format", default="GT", help="format tag to compare (default GT)")
    (options, args)=parser.parse_args()
    vcfilename=args[0]
    #vcfilename='/Users/indapa/software/Pgmsnp/PythonNotebook/child5x.nrs.sites.calledWith20x_bam.child5x.nrs.sites.calledWith5x_bam.combineVariants.vcf'
    
    basename=os.path.splitext(vcfilename)[0]

    vcfobj=VcfFile(vcfilename)
    vcfh=open(vcfilename,'r')

    vcfobj.parseMetaAndHeaderLines(vcfh)
    header=vcfobj.returnHeader() +"\n"

    samples=vcfobj.getSampleList()
    print "\t".join(samples)
    for vrec in vcfobj.yieldVcfRecordwithGenotypes(vcfh):
        
        vrec_ziptuple=vrec.zipGenotypes(samples)
        outputs=[]
        for (sample, geno_obj) in vrec_ziptuple:
            outputs.append( "\t".join( [geno_obj.getFormatVal(options.format) ] ) )
        print "\t".join(outputs)


if __name__ == "__main__":
    main()
