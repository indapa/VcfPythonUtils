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



def computeCalledGenotypes (gtm):
    """ given a genotype matrix of counts of genotypes in the 4 possible classes, sum the total called genotypes """

    total_calls = gtm[0,0] + gtm[0,1] + gtm[0,2] + gtm[1,0] + gtm[1,1] + gtm[1,2] + gtm[2,0] + gtm[2,1] + gtm[2,2]
    return total_calls

def computeRowMarginalAltCount (gtm):
    """ given a genotype matrix of counts of genotypes in the 4 possible classes,  compute the row marginal alt allele count  """
    het_calls = gtm[1,0] + gtm[1,1] + gtm[1,2] + gtm[1,3]
    homoz_alt_calls = gtm[2,0] + gtm[2,1] + gtm[2,2]  + gtm[2,3]
    alt_count = het_calls + (2*homoz_alt_calls)
    return alt_count

def computeColMarginalAltCount (gtm):
    """ given a genotype matrix of counts of genotypes in the 4 possible classes,  compute the column marginal alt allele count  """
    het_calls = gtm[0,1] + gtm[1,1] + gtm[2,1] + gtm[3,1]
    homoz_alt_count = gtm[0,2] + gtm[1,2] + gtm[2,2] + gtm[3,2]

    alt_count = het_calls + (2*homoz_alt_count)
    return alt_count

def computeNRD(gtm):
    """ compute the Non-reference discrepancy rate: http://www.broadinstitute.org/gsa/wiki/index.php/File:GenotypeConcordanceGenotypeErrorRate.png  """
    """ ignores concordant calls that are  ref/ref; conditions on calls being made by both callsets """


    discordant_call_count = gtm[0,1]+gtm[0,2]+gtm[1,0]+gtm[1,2]+gtm[2,0]+gtm[2,1]
    total_count = gtm[0,1]+gtm[0,2]+gtm[1,0]+gtm[1,1]+ gtm[1,2]+gtm[2,0]+gtm[2,1] +gtm[2,2]

    if total_count==0: return 'NA'

    nrd= round( float(discordant_call_count)/float(total_count), 3)
    return nrd

def computeNRD_class (gtm):
    """ compute the discrepancy in each genotype class (homoz_ref, het, homoz_nonref) """
    total_homoz_ref= gtm[0,0] + gtm[1,0] + gtm[2,0]
    discord_homoz_ref= gtm[1,0] + gtm[2,0]

    total_het= gtm[0,1] + gtm[1,1] + gtm[2,1]
    discord_het= gtm[0,1] + gtm[2,1]

    total_homoz_nonref= gtm[0,2] + gtm[1,2] + gtm[2,2]
    discord_homoz_nonref = gtm[0,2] + gtm[1,2]

    if total_homoz_ref == 0:
        nrd_homoz_ref='NA'
    else:
        nrd_homoz_ref= round ( float(discord_homoz_ref)/float(total_homoz_ref), 3 )

    if total_het == 0:
        nrd_het = 'NA'
    else:
        nrd_het = round( float(discord_het)/float(total_het), 3 )

    if total_homoz_nonref == 0:
        nrd_homoz_ref = 'NA'
    else:
        nrd_homoz_ref= round ( float(discord_homoz_nonref)/float(total_homoz_nonref), 3 )

    return (nrd_homoz_ref, nrd_het, nrd_homoz_ref,
            (total_homoz_ref, discord_homoz_ref),
            (total_het, discord_het),
            (total_homoz_nonref, discord_homoz_nonref) )


def computeNRS(gtm):
    """ compute non-reference sensitivity: http://www.broadinstitute.org/gsa/wiki/images/0/07/GenotypeConcordanceVarSens.png """
    """ number variant genotypes in evalutationset / number of variant genotypes in comparison truthset """

    variant_count_evaluation= gtm[1,1]+ gtm[1,2]+ gtm[2,1]+ gtm[2,2]

    variant_count_comparison= gtm[0,1]+gtm[0,2]+gtm[1,1]+gtm[1,2]+gtm[2,1]+gtm[2,2]+gtm[3,1]+gtm[3,2]

    if variant_count_comparison ==0:
        return 'NA'
    nrs= round ( float(variant_count_evaluation)/float(variant_count_comparison), 3 )
    return nrs


def binGenoytpeQualities( gq_calibrations):
    """ given a calibrationList [ (gq, 0|1), ... ] bin the data"""
    bins = np.linspace(0, 100, 11) #11 bins between 0 and 100 (0,10,20,30,40,50,60,70,80,90,100)
    L =[ [ ] for i in range(0, len(bins)) ] # L is a list of lists [ [] ,,, [] ]
    outfh=open('gq.calibration.txt', 'w')

    for (gq, truth) in gq_calibrations:

        for i in range(0, len(bins)-1 ):
            if int(gq) >= bins[i] and int(gq) < bins[i+1]:
                L[i].append(  truth )
                #print str(gq), truth
                continue
    for i in range(0, len(L) ):
        total=len(L[i])
        totalinCorrect=0
        for flag in L[i]:
            if flag == 0: totalinCorrect+=1

        j=bins[i]
        outstring = "\t".join([str(j), str(totalinCorrect), str(total)])
        outfh.write(outstring+"\n")

    outfh.close()

    

def main():
    """compare the genotypes in second vcf file to the ones in the first vcf file
       the two VCFs should have the same sites but not necessarily the same samples
       output is the concordance rate of the genotypes of the samples in the first VCF to the genotypes of the samples in the second
       
       if there are no common samples between the two VCFs no output is written!
    """
    description="\ncompare the genotypes in second vcf file to the ones in the first vcf file\nThe two VCFs should have the same sites but not necessarily the same samples"
    usage = "usage: %prog file1.vcf file2.vcf\n"+description

    parser = OptionParser(usage)
    parser.add_option("--onlyImputed",  action="store_true", dest="imputedonly", default=False, help="set option if you want to compare *only imputed* genotypes in the first file  to the genotypes in the second")
    parser.add_option("--ignoreImputed",  action="store_true", dest="ignoreimputed", default=False, help="set option if you want to *ignore* imputed genotypes in the first file  to the genotypes in the second")
    parser.add_option("--gprob", type="float", dest="gprob", default=0.0, help="for comparison to genotypes in second VCF, *imputed* genotype probabilties in first VCF must be at least gprob and have the format tag GPROB")
    parser.add_option("--r2", type="float", dest="rsquare", default=0.0, help="for comparison of *imputed*  genotypes in first VCF to those  in second VCF, SNP r2 of the site in first VCF must be at least rsquare and have the  info tag R2")


    (options, args)=parser.parse_args()
    (vcf_fname1, vcf_fname2) =  args[0:2]
   
    
    vcf_fh1= open(vcf_fname1, 'r')
    vcf_fh2= open(vcf_fname2, 'r')
 
    logfh= open('mismatch.log', 'w')
    sitefh = open('site.nrd.nrs.txt', 'w')
    site_header= "\t".join( [ 'chrom', 'pos', 'ref', 'alt', 'rsquare', 'total', 'called', 'totalChrom', 'vcf1_ac', 'vcf1_nocall', 'vcf2_ac', 'vcf2_nocall', 'vcf2_maf', 'nrs', 'nrd', 'nrd_homoz_ref', 'nrd_het', 'nrd_homoz_nonref'])
    sitefh.write(site_header+'\n')

    vcf1_samples=get_vcfsamples(vcf_fh1)
    vcf2_samples=get_vcfsamples(vcf_fh2)


    sitetotalAA=0
    sitetotalAA_discord=0

    sitetotalAB=0
    sitetotalAB_discord=0

    sitetotalBB=0
    sitetotalBB_discord=0


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
   
    genotype_quality_list=[] # [ (gq, 0|1) .... ]


    #reset the filehandle positions
    vcf_fh1.seek(0)
    vcf_fh2.seek(0)

    
    while 1:
        if '#CHROM' in vcf_fh1.readline(): break

    while 1:
        if '#CHROM' in vcf_fh2.readline(): break

   
    if options.imputedonly == True:
        sys.stderr.write("comparing only imputed genotypes in first VCF to genotypes in the second\n")
    else:
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
        vcf1_infostr = vcf1_data[7]
        #print vcf1_formatstr
        
        #print vcf1_data[0], vcf1_data[1],  vcf1_data[3:5], vcf2_data[3:5]

        #the genotype comparison matrix for the VCF site
        #reset everytime a new site is analyzed
        site_discordance=np.matrix( [ [ 0,0,0,0 ], [ 0,0,0,0 ], [ 0,0,0,0 ], [ 0,0,0,0 ] ] )

        #print site_discordance

        r2value='NA'
        if 'R2' in vcf1_infostr:
            infofields=vcf1_infostr.split(';')
            #print infofields
            if 'R2' in infofields[3]:
                (r2,r2value)=infofields[3].split('=')
                value=float(r2value)
                if float(r2value) <= options.rsquare:
                    #print options.rsquare, value
                    continue

        if vcf1_data[3:5] != vcf2_data[3:5]:
            sys.stderr.write("ref/alt alleles don't match between vcfs!\n")
            outstr = "\t".join( [vcf1_data[0], vcf1_data[1], vcf1_data[3], vcf1_data[4], vcf2_data[0], vcf2_data[1], vcf2_data[3], vcf2_data[4] ] )
            logfh.write(outstr+"\n")
            continue

        if vcf1_data[0:2] != vcf2_data[0:2]:
            sys.stderr.write("chrom/position doesn't match!")
            print vcf1_data[0], vcf1_data[1],vcf2_data[0], vcf2_data[1]
            exit(1)
        

        #ziptuple=zip[ (sample, genotypefield), .... ]

        vcf1_ziptuple=zip(vcf1_samples, vcf1_data[9::] )
        vcf2_ziptuple=zip(vcf2_samples, vcf2_data[9::] )

        #filetr the ziptuple to contain only those samples common to both!
        filtered_vcf1= [x for x in vcf1_ziptuple if x[0] in common_samples]
        filtered__vcf2=  [x for x in vcf2_ziptuple if x[0] in common_samples]

        maf_vcf2 = calMaf(filtered__vcf2) #get the MAF for the site in the second VCF
        #print len(filtered_vcf1), len(filtered__vcf2)
       
        #collect the compariosn results
        if (options.imputedonly == True):
            comparison_results= compare_imputed_genotypes(filtered_vcf1, filtered__vcf2, vcf1_formatstr, options.gprob)
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
            
        elif (options.ignoreimputed == True):
            comparison_results= compare_nonimputed_genotypes(filtered_vcf1, filtered__vcf2, vcf1_formatstr, options.gprob)
        else:
            comparison_results = compare_genotypes(filtered_vcf1, filtered__vcf2)
            gq_calibrations = gq_calibration(filtered_vcf1, filtered__vcf2, vcf1_formatstr)
            for t in gq_calibrations:
                genotype_quality_list.append(t)
            #print comparison_results

        for (g1, g2, sample) in comparison_results: #iterate thru the compariosn results
            discordance_dict[sample][g1,g2]+=1      #incrment the per-sample counts of genotype compariosns results
            site_discordance[g1,g2]+=1              #incrment the per-site counts of genotype compariosns results

        #print site_discordance
        #collect site specific nrd, nrs, called genotypes, alt count, etc for this current site
        site_nrs=computeNRS( site_discordance )
        site_nrd=computeNRD ( site_discordance )
        (site_nrd_homoz_ref, site_nrd_het, site_nrd_homoz_nonref, aa, ab, bb )=computeNRD_class( site_discordance)

        sitetotalAA+=aa[0]
        sitetotalAA_discord+=aa[1]

        sitetotalAB+=ab[0]
        sitetotalAB_discord+=ab[1]

        sitetotalBB+=bb[0]
        sitetotalBB_discord+=bb[1]


        #if site_nrd == 'NA': print site_discordance
        vcf1_nocalls = site_discordance[3,0] + site_discordance[3,1] + site_discordance[3,2] + site_discordance[3,3]
        vcf2_nocalls = site_discordance[0,3] + site_discordance[1,3] + site_discordance[2,3] + site_discordance[3,3]


        totalCalledGenotypes = computeCalledGenotypes(site_discordance)
        totalChroms=2*totalCalledGenotypes
        
        vcf1_altcount=computeRowMarginalAltCount(site_discordance)
        vcf2_altcount =computeColMarginalAltCount(site_discordance)

        #if vcf1_altcount == 1: print "ac ==1\n", site_discordance
        #if site_nrd == 1.0: print "nrd==1\n", site_discordance

        outstr= "\t".join( [ vcf1_data[0], vcf1_data[1], vcf1_data[3],
                             vcf1_data[4], str(r2value), str( len(filtered_vcf1) ),
                             str(totalCalledGenotypes), str(totalChroms), str(vcf1_altcount),
                             str(vcf1_nocalls), str(vcf2_altcount), str(vcf2_nocalls) ])

        sitefh.write(outstr +"\t")
        sitefh.write('%.3f'%maf_vcf2+'\t')
        sitefh.write( str(site_nrs)+'\t')
        sitefh.write( str(site_nrd)+'\t')
        sitefh.write( str(site_nrd_homoz_ref) + '\t')
        sitefh.write( str(site_nrd_het) + '\t')
        sitefh.write( str(site_nrd_homoz_nonref) + '\n')
       

    #broke out of the while loop - finished comparing the two vcf files
    #now collect and write per sample nrd and nrs results
    outfh=open('sample.nrd.nrs.txt', 'w')
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
        outstring = "\t".join( [ sample, str(nrd), str(nrc),  str ( np.sum( discordance_dict[sample] ) ) , str(eval_nocalls) , str(comparison_nocalls ) ] )
        outfh.write(outstring+"\n")

    outfh.close()

    #sys.stderr.write("writing accuracy v. posterior prob calibration...\n")
    #calibrateBins(calibration_matrix)
    #calibmatrix = binGenoytpeQualities(genotype_quality_list)

    outstr="\t".join( [str(sitetotalAA), str(sitetotalAA_discord),
                      str(sitetotalAB), str(sitetotalAB_discord),
                      str(sitetotalBB), str(sitetotalBB_discord)] )
    print outstr

if __name__ == "__main__":
    main()