#!/usr/bin/python
import sys
import os
import string
import re
from optparse import OptionParser
from vcfIO import *
import numpy as np


def calibrateBins (matrix):
    outfh=open('calibration.txt', 'w')
    for i in range(0, len(matrix) ):
        total=len(matrix[i])
        totalCorrect=0
        for (prob, flag) in matrix[i]:
            if flag == 1: totalCorrect+=1
        j=float(i)/10
        outstring = "\t".join([str(j), str(totalCorrect), str(total)])
        outfh.write(outstring+"\n")

    outfh.close()

def computeNRD(gtm):
    """ compute the Non-reference discrepancy rate: http://www.broadinstitute.org/gsa/wiki/index.php/File:GenotypeConcordanceGenotypeErrorRate.png  """
    """ ignores concordant calls that are  ref/ref; conditions on calls being made by both callsets """


    discordant_call_count = gtm[0,1]+gtm[0,2]+gtm[1,0]+gtm[1,2]+gtm[2,0]+gtm[2,1]
    total_count = gtm[0,1]+gtm[0,2]+gtm[1,0]+gtm[1,1]+ gtm[1,2]+gtm[2,0]+gtm[2,1] +gtm[2,2]

    if total_count==0: return 'NA'

    nrd= float(discordant_call_count)/float(total_count)
    return nrd


def computeNRS(gtm):
    """ compute non-reference sensitivity: http://www.broadinstitute.org/gsa/wiki/images/0/07/GenotypeConcordanceVarSens.png """
    """ number variant genotypes in evalutationset / number of variant genotypes in comparison truthset """

    variant_count_evaluation= gtm[1,1]+ gtm[1,2]+ gtm[2,1]+ gtm[2,2]

    variant_count_comparison= gtm[0,1]+gtm[0,2]+gtm[1,1]+gtm[1,2]+gtm[2,1]+gtm[2,2]+gtm[3,1]+gtm[3,2]

    if variant_count_comparison ==0:
        return 'NA'
    nrs= float(variant_count_evaluation)/float(variant_count_comparison)
    return nrs


def binCalibrations( gprobs_calibrations, calibration_matrix):
    """ given a calibrationList [ (gprob, 0|1), ... ] and a 2d matrix calibrationMatrix append the tuple to the proper bin (row) in the matrix """

    
    return calibration_matrix

def main():
    """compare the genotypes in second vcf file to the ones in the first vcf file
       the two VCFs should have the same sites but not necessarily the same samples
       output is the concordance rate of the genotypes of the samples in the first VCF to the genotypes of the samples in the second
       
       if there are no common samples between the two VCFs no output is written!
    """
    description="\ncompare the genotypes in second vcf file to the ones in the first vcf file\nThe two VCFs should have the same sites but not necessarily the same samples"
    usage = "usage: %prog file1.vcf file2.vcf\n"+description

    parser = OptionParser(usage)
    parser.add_option("--compareImputed",  action="store_true", dest="imputedonly", default=False, help="set option if you want to compare only imputed genotypes in the first file  to the genotypes in the second")
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

    
    calibration_matrix= []
    for i in range(0,10):
        calibration_matrix.append( [] )
   



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
     
        vcf1_formatstr= vcf1_data[8]

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
        if (options.imputedonly == True):
            comparison_results= compare_imputed_genotypes(filtered_vcf1, filtered__vcf2, vcf1_formatstr)
            gprobs_calibrations = posterior_imputed_gprob_calibration(filtered_vcf1, filtered__vcf2, vcf1_formatstr)
            for x in gprobs_calibrations:
                (gprob, comparison_result) =x
                gprob =float(gprob)
                
                if gprob <= .10:
                    calibration_matrix[0].append(x)
                elif gprob <= .20:
                    calibration_matrix[1].append(x)
                elif gprob <= .30:
                    calibraton_matrix[2].append(x)
                elif gprob <= .40:
                    calibration_matrix[3].append(x)
                elif gprob <= .50:
                    calibration_matrix[4].append(x)
                elif gprob <= .60:
                    calibration_matrix[5].append(x)
                elif gprob <= .70:
                    calibration_matrix[6].append(x)
                elif gprob <= .80:
                    calibration_matrix[7].append(x)
                elif gprob <= .90:
                    calibration_matrix[8].append(x)
                elif gprob <= 1.0:
                    calibration_matrix[9].append(x)
                else:
                    pass
            


        else:
            comparison_results = compare_genotypes(filtered_vcf1, filtered__vcf2)

        for (g1, g2, sample) in comparison_results:
            discordance_dict[sample][g1,g2]+=1



    outfh=open('nrd.nrs.txt', 'w')
    headerstring = "\t".join( [ "sample","NRD", "NRS", "totalGenotypes", "noCallsEval", "noCallsComparison"] )
    outfh.write(headerstring+"\n")
    for sample in discordance_dict.keys():
        #print sample
        #print discordance_dict[sample]
        #print "total gneotypes: ", np.sum( discordance_dict[sample] )

        nrc=computeNRS( discordance_dict[sample] )
        nrd=computeNRD(  discordance_dict[sample]  )
        gtm= discordance_dict[sample]
        eval_nocalls = gtm[3,0] + gtm[3,1] + gtm[3,2] + gtm[3,3]
        comparison_nocalls = gtm[0,3] + gtm[1,3] + gtm[2,3] + gtm[3,3]
        
        outstring = "\t".join(sample, str(nrd), str(nrc),  str ( np.sum( discordance_dict[sample] ) ) , str(eval_nocalls) , str(comparison_nocalls) )
        outfh.write(outstring+"\n")

    outfh.close()

    sys.stderr.write("writing accuracy v. posterior prob calibration...\n")
    calibrateBins(calibration_matrix)

if __name__ == "__main__":
    main()