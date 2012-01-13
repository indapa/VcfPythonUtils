#!/usr/bin/python
import sys
import os
import string
import re
from optparse import OptionParser

from VcfFile import *

""" look for shared variants in a list of affecteds given an inheritance model """


def main():
    usage = "usage: %prog [options] file.vcf"
    parser = OptionParser(usage)
    parser.add_option("--model", type="string", dest="model", default = "dominant", help=" inheritance model (dominant recessive) ")
    parser.add_option("--affected", type="string", dest="affected", help="sample name of affecteds (one per line)")
    parser.add_option("--filter", type="string", dest="filter", help="analyze only those  records matching filter (default is None)", default=None)
    (options, args)=parser.parse_args()
    if options.affected==None:
        sys.stderr.write("please provide a value to --affected parameter!\n")
        exit(1)

    affecteds=[] # list of affected samples
    affectedfh=open(options.affected, 'r')
    for line in affectedfh:
        affecteds.append(line.strip() )

    vcfilename=args[0]
    vcfh=open(vcfilename,'r')

    #instantiate a VcfFile object
    vcfobj=VcfFile(vcfilename)
    vcfobj.parseHeaderLine(vcfh)
    samplelist=vcfobj.getSampleList()
    print "samplelist", samplelist
    vcfh.seek(0)
    #parse its metainfo lines (ones that begin with ##)
    vcfobj.parseMetaLines(vcfh)

    for vrec in vcfobj.yieldVcfRecordwithGenotypes(vcfh ):
        
        affected_genotypes=[] #list of tuples (sample, genotypeobj) with samples that are affected
        if vrec.getFilter() != options.filter and options.filter != None : continue
        genotypes = vrec.getGenotypesAlleles()
        genotype_tuple= vrec.zipGenotypes(samplelist)
        for (sample, genotype) in genotype_tuple:
            if sample in affecteds:
                affected_genotypes.append( ( sample, genotype.isSegregating() )  )
        #print affected_genotypes
        shared_affected_segregating = filter( lambda x, segregating=True: segregating in x, affected_genotypes)
        #now if all affects are segregating for the site, its a candidate
        #print len(shared_affected_segregating), len(affecteds), shared_affected_segregating, affected_genotypes
        if len(shared_affected_segregating) == len(affecteds):
            print shared_affected_segregating
if __name__ == "__main__":
    main()
