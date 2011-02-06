#!/usr/bin/python
import sys
import os
import string
import re
from optparse import OptionParser
from vcfIO import *

def main():
    """compare the genotypes in second vcf file to the ones in the first vcf file
       the two VCFs should have the same sites but not necessarily the same samples
       output is the concordance rate of the genotypes of the samples in the first VCF to the genotypes of the samples in the second
       
       if there are no common samples between the two VCFs no output is written!
    """
    description="\ncompare the genotypes in second vcf file to the ones in the first vcf file\nThe two VCFs should have the same sites but not necessarily the same samples"
    usage = "usage: %prog file1.vcf file2.vcf\n"+description

    parser = OptionParser(usage)
    
    (options, args)=parser.parse_args()
    (vcf_fname1, vcf_fname2) =  args[0:2]
   
    
    vcf_fh1= open(vcf_fname1, 'r')
    vcf_fh2= open(vcf_fname2, 'r')
 
    

    vcf1_samples=get_vcfsamples(vcf_fh1)
    vcf2_samples=get_vcfsamples(vcf_fh2)


    common_samples= []
    for s in vcf2_samples:
        if s in vcf1_samples: common_samples.append(s)

    if len(common_samples) == 0:
        sys.stderr.write("vcf files have no common samples!\n")
        exit(1)

    #key samplename value discordance count of genotypes between the two vcf files
    discordance_dict={}
    for s in common_samples:
        discordance_dict[s]=0

    
    #reset the filehandle positions
    vcf_fh1.seek(0)
    vcf_fh2.seek(0)

    
    while 1:
        if '#CHROM' in vcf_fh1.readline(): break

    while 1:
        if '#CHROM' in vcf_fh2.readline(): break

   
    sys.stderr.write("comparing genotypes between vcf files...\n")


    while 1:
       
        vcf1_line= vcf_fh1.readline()
        vcf2_line =  vcf_fh2.readline()

        if vcf1_line == '' or vcf2_line == '':
            #sys.stderr.write("unphased vcf has reached EOF!\n")
            break
        
        vcf1_data=split_vcfdataline(vcf1_line)
        vcf2_data=split_vcfdataline(vcf2_line)
     
        if vcf1_data[0:2] != vcf2_data[0:2]:
            sys.stderr.write("chrom/position doesn't match!")
            print vcf1_data[0:2], vcf2_data[0:2]
            exit(1)
        
        #ziptuple=zip[ (sample, genotypefield), .... ]

        vcf1_ziptuple=zip(vcf1_samples, vcf1_data[9::] )
        vcf2_ziptuple=zip(vcf2_samples, vcf2_data[9::] )

        #filetr the ziptuple to contain only those samples common to both!
        filtered_vcf1= [x for x in vcf1_ziptuple if x[0] in common_samples]
        filtered__vcf2=  [x for x in vcf2_ziptuple if x[0] in common_samples]

        
        #collect the unmatched genotypes at a locus, and if there are any, print to STDOUT
        unmatched_genotypes = compare_genotypes(filtered_vcf1, filtered__vcf2)
        if len(unmatched_genotypes) > 0:
            #print vcf1_data[0], vcf1_data[1], "unmatched genotypes: ", unmatched_genotypes
            for (g1,g2,sample) in unmatched_genotypes:
                discordance_dict[sample]+=1
            #print "\n"
    for s in discordance_dict.keys():
        print s, discordance_dict[s]

if __name__ == "__main__":
    main()