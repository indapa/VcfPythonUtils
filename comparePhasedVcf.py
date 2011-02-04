#!/usr/bin/python
import sys
import os
import string
import re
from optparse import OptionParser
from vcfIO import *

def main():
    """ given the unphased and phased version of VCF file, check for genotype inconistencies between phased genotypes and unphased ones """
    usage = "usage: %prog [options] phased.vcf"

    parser = OptionParser(usage)
    parser.add_option("--unphasedvcf", type="string", dest="unphased", help="unphased vcf ")

    (options, args)=parser.parse_args()
    phasedvcf_file=args[0]
   
    
    unphased_fh=open(options.unphased, 'r')
    

    phased_fh=open(phasedvcf_file, 'r')
    
    

    phased_samples=get_vcfsamples(phased_fh)
    unphased_samples=get_vcfsamples(unphased_fh)

    #print phased_samples
    #print unphased_samples

    #compare the sample lists to makesure the samples are the same
    if compareSampleLists(phased_samples, unphased_samples) == 0:
        sys.stderr.write("vcf files have different samples!")
        exit(1)

    #the phased adn unphased vcf should have teh same number of loci!
    phased_fh.seek(0)
    unphased_fh.seek(0)
    

    #for phasedline in phased_fh:
    while 1:
       #print line.strip()
        unphasedline=unphased_fh.readline()
        phasedline = phased_fh.readline()

        if unphasedline == '' or phasedline == '':
            sys.stderr.write("unphased vcf has reached EOF!")
            exit(1)
        
        
        
        if '#' in phasedline or '#' in unphasedline: continue


        #print phasedline
        #print unphasedline

        phased_data=split_vcfdataline(phasedline)
        unphased_data=split_vcfdataline(unphasedline)

        #check to see if the chrom/pos match between the phased adn unphased
        if phased_data[0:2] != unphased_data[0:2]:
            sys.stderr.write("chrom/position doesn't match!")
            print phased_data[0::2], unphased_data[0:2]
            exit(1)
        
        #ziptuple=zip[ (sample, genotypefield), .... ]
        phased_ziptuple=zip(phased_samples, phased_data[9::] )
        unphased_ziptuple=zip(unphased_samples, unphased_data[9::] )

        #collect the unmatched genotypes at a locus, and if there are any, print to STDOUT
        unmatched_genotypes =compare_phased_to_unphased(phased_ziptuple, unphased_ziptuple)
        if len(unmatched_genotypes) > 0:
            print phased_data[0], phased_data[1], "unmatched genotypes: ", unmatched_genotypes
        

if __name__ == "__main__":
    main()