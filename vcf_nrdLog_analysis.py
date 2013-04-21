#!/usr/bin/env python
from itertools import *
from VcfFile import *
import numpy as np
import re, os
from optparse import OptionParser
from common import grouper, typeofGenotype, melt_lol

""" do some data massaging to figure out patterns in sites that have
discrepant genotypes. Input file is the vcf logfile of sites that contribute
to NRD in the input file to variantEval """

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
    nrdallfh=open(basename+".allgenos.nrd.txt", 'w')
    nrdtwofh=open(basename+".twogenos.nrd.txt", 'w')
    nrdonefh=open(basename+".onegenos.nrd.txt", 'w')
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
        """ Since I'm testing against trio, NRD count can be 1 2 or 3
            We keep track of the nrd count and print those records to the appropriate file:
            nrdallfh, nrdtwofh, nrdonefh  """
        nrd_count=0
        for (compare, eval) in grouper(2,vrec_ziptuple):
            (comp_allele1, comp_allele2)=compare[1].getAlleles()
            (eval_allele1, eval_allele2)=eval[1].getAlleles()
            eval_alleletype=typeofGenotype(eval_allele1, eval_allele2)
            comp_alleletype=typeofGenotype(comp_allele1, comp_allele2)
            if eval_alleletype == comp_alleletype:
                continue
            outputline.append(   [ eval[0], str(eval_alleletype), compare[0], str(comp_alleletype) ] )
            nrd_count+=1
        
        output= "\t".join (  melt_lol(outputline) )
        """ depending on the nrd count, print the records to appropirate file(s) """
        if nrd_count == 3:
            nrdallfh.write(output+"\n")
        if nrd_count== 2:
            nrdtwofh.write(output+"\n")
        if nrd_count==1:
            nrdonefh.write(output+"\n")
            
            
       
            
        
        

if __name__ == "__main__":
    main()