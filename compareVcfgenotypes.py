#!/usr/bin/python
import sys
import os
import string
import re
from optparse import OptionParser
from vcfIO import *
import numpy as np


def computeNRD(gtm):
    """ compute the Non-reference discrepancy rate: http://www.broadinstitute.org/gsa/wiki/index.php/File:GenotypeConcordanceGenotypeErrorRate.png  """
    """ ignores concordant calls that are  ref/ref; conditions on calls being made by both callsets """


    discordant_call_count = gtm[0,1]+gtm[0,2]+gtm[1,0]+gtm[1,2]+gtm[2,0]+gtm[2,1]
    total_count = gtm[0,1]+gtm[0,2]+gtm[1,0]+gtm[1,1]+ gtm[1,2]+gtm[2,0]+gtm[2,1] +gtm[2,2]
    nrd= float(discordant_call_count)/float(total_count)
    return nrd


def computeNRS(gtm):
    """ compute non-reference sensitivity: http://www.broadinstitute.org/gsa/wiki/images/0/07/GenotypeConcordanceVarSens.png """
    """ number variant genotypes in evalutationset / number of variant genotypes in comparison truthset """

    variant_count_evaluation= gtm[1,1]+ gtm[1,2]+ gtm[2,1]+ gtm[2,2]

    variant_count_comparison= gtm[0,1]+gtm[0,2]+gtm[1,1]+gtm[1,2]+gtm[2,1]+gtm[2,2]+gtm[3,1]+gtm[3,2]

    nrs= float(variant_count_evaluation)/float(variant_count_comparison)
    return nrs



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
        discordance_dict[s]=np.matrix( [ [ 0,0,0,0 ], [ 0,0,0,0 ], [ 0,0,0,0 ], [ 0,0,0,0 ] ] )

    
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

        
        
        #collect the compariosn results
        comparison_results = compare_genotypes(filtered_vcf1, filtered__vcf2)

        for (g1, g2, sample) in comparison_results:
            discordance_dict[sample][g1,g2]+=1
    print "sample","NRD", "NRS", "totalGenotypes", "noCallsEval", "noCallsComparison"
    for sample in discordance_dict.keys():
        #print sample
        #print discordance_dict[sample]
        #print "total gneotypes: ", np.sum( discordance_dict[sample] )

        nrc=computeNRS( discordance_dict[sample] )
        nrd=computeNRD(  discordance_dict[sample]  )
        missing_genotypes = discordance_dict[sample][3,0]
        missing_genotypes_comparison= discordance_dict[sample][0,3]

        print sample,nrd, nrc,  np.sum( discordance_dict[sample] ) , np.sum( missing_genotypes ), np.sum( missing_genotypes_comparison )
        #print "=="
if __name__ == "__main__":
    main()