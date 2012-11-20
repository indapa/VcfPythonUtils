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
    parser.add_option("--model", type="string", dest="model", default = "dominant", help=" inheritance model (dominant/recessive), default is dominant ")
    parser.add_option("--affected", type="string", dest="affected", help="sample name of affecteds (one per line)")
    parser.add_option("--unaffected", type="string", dest="unaffected", help="sample name of unaffecteds (one per line)")
    parser.add_option("--filter", type="string", dest="filter", help="analyze only those  records matching filter (default is None)", default=None)
    (options, args)=parser.parse_args()
    if options.affected==None:
        sys.stderr.write("please provide a value to --affected parameter!\n")
        exit(1)

    affecteds=[] # list of affected samples
    unaffecteds=[] # list of unaffected samples

    affectedfh=open(options.affected, 'r')
    for line in affectedfh:
        affecteds.append(line.strip() )

    if options.unaffected != None:
        unaffectedfh=open (options.unaffected, 'r')
        for line in unaffectedfh:
            unaffecteds.append ( line.strip() )

    #check if any overlapping samples between unaffected and affected
    if len( list( set(unaffecteds).intersection( set(affecteds) ) )  ) != 0:
        sys.stderr.write("check list of affected and unaffecteds for overlapping samples!\n")
        exit(1)

    #    sys.stderr.write("check list of affected and unaffected for overlapping samples!\n")
    #    exit(1)


    vcfilename=args[0]
    vcfh=open(vcfilename,'r')

    #instantiate a VcfFile object
    vcfobj=VcfFile(vcfilename)
    #print the ## header lines
    vcfobj.parseMetaLines(vcfh)
    vcfobj.printMetaLines()
    vcfh.seek(0)

    #part the #Chrom header line
    vcfobj.parseHeaderLine(vcfh)
    vcfobj.printHeaderLine()
    samplelist=vcfobj.getSampleList()
    
    
    for vrec in vcfobj.yieldVcfRecordwithGenotypes(vcfh ):
        
        affected_genotypes=[] #list of tuples (sample, VcfGenotype object) with samples that are affected
        unaffected_genotypes=[] # list of tuples (sample, VcfGenotype object) with samples that are unaffected

        if vrec.getFilter() != options.filter and options.filter != None : continue
        genotypes = vrec.getGenotypesAlleles()
        genotype_tuple= vrec.zipGenotypes(samplelist) # get a list of tuples [ (sample, VcfGenotype object) ... ]
        for (sample, genotype) in genotype_tuple: #iterate thru and see if they are in affected or unaffected list
            if options.model == 'dominant':
                if sample in affecteds:  # if so ...
                    affected_genotypes.append( ( sample, genotype.toString(),  genotype.isSegregating() )  ) # are they segregating for a non-ref allele?
                if sample in unaffecteds:
                    unaffected_genotypes.append( (sample,  genotype.toString(),  genotype.isSegregating() ) ) # are they segregating for a non-ref allele?
            elif options.model == 'recessive':
                if sample in affecteds:
                    affected_genotypes.append( ( sample, genotype.toString(),  genotype.isNonRefHomz() )  ) # are they segregating for a non-ref homoz?
                if sample in unaffecteds:
                    unaffected_genotypes.append( (sample,  genotype.toString(),  genotype.isNonRefHomz() ) ) # are they segregating for a non-ref non-refhomoz?
            else:
                sys.stderr.write(options.model + " not supported for genotype discrete filtering ...\n")


        if options.model == 'dominant':
            #filter the collected samples to see if they are all have segregating genotypes
            shared_affected_segregating = filter( lambda x, segregating=True: segregating in x, affected_genotypes)
            shared_unaffected_segregating = filter ( lambda x, segregating=False: segregating in x, unaffected_genotypes)
        
            #now if all affects are segregating for the site
            # and all the un-affecteds are *not* segregating for the site
            # it is a candidate
            if len(shared_affected_segregating) == len(affecteds):
                if  len(shared_unaffected_segregating) == len(unaffecteds):
                    print vrec.toStringwithGenotypes()
        elif options.model == 'recessive':
            shared_affected_homoz = filter( lambda x, homoz=True: homoz in x, affected_genotypes)
            shared_unaffected_homoz = filter ( lambda x, homoz=False: homoz in x, unaffected_genotypes)

            #now if all affects are homoz for the site
            # and all the un-affecteds are *not* homoz for the site
            # it is a candidate
            if len(shared_affected_homoz) == len(affecteds):
                if  len(shared_unaffected_homoz) == len(unaffecteds):
                    print vrec.toStringwithGenotypes()
        else:
            sys.stderr.write(options.model + " not supported for genotype discrete filtering ...\n")

        #print "\n"
if __name__ == "__main__":
    main()
