#!/usr/bin/python
import sys
import os
import string
import re
from optparse import OptionParser
from vcfIO import *
import numpy as np

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



def main():
    usage = "usage: %prog [options] file1.vcf file2.vcf\n\nCompare heterzygote genotype discrepancies in the fist VCF to the second VCF"
    parser = OptionParser(usage)
    #parser.add_option("--name", type="string|int|boolean", dest="name", help="help string")

    (options, args)=parser.parse_args()

    (vcf_fname1, vcf_fname2) =  args[0:2]
   
    
    vcf_fh1= open(vcf_fname1, 'r')
    vcf_fh2= open(vcf_fname2, 'r')
    
    vcf1_samples=get_vcfsamples(vcf_fh1)
    vcf2_samples=get_vcfsamples(vcf_fh2)
    common_samples=[]
    
    for s in vcf2_samples:
        if s in vcf1_samples: common_samples.append(s)

    if len(common_samples) == 0:
        sys.stderr.write("vcf files have no common samples!\n")
        exit(1)
        
     #reset the filehandle positions
    vcf_fh1.seek(0)
    vcf_fh2.seek(0)

    
    while 1:
        if '#CHROM' in vcf_fh1.readline(): break

    while 1:
        if '#CHROM' in vcf_fh2.readline(): break
        
    
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
        

        #the genotype comparison matrix for the VCF site
        #reset everytime a new site is analyzed
        site_discordance=np.matrix( [ [ 0,0,0,0 ], [ 0,0,0,0 ], [ 0,0,0,0 ], [ 0,0,0,0 ] ] )
        

        if vcf1_data[3:5] != vcf2_data[3:5]:
            sys.stderr.write("ref/alt alleles don't match between vcfs!\n")
            outstr = "\t".join( [vcf1_data[0], vcf1_data[1], vcf1_data[3], vcf1_data[4], vcf2_data[0], vcf2_data[1], vcf2_data[3], vcf2_data[4] ] )
            logfh.write(outstr+"\n")
            continue

        if vcf1_data[0:2] != vcf2_data[0:2]:
            sys.stderr.write("chrom/position doesn't match!")
            print vcf1_data[0], vcf1_data[1],vcf2_data[0], vcf2_data[1]
            exit(1)
        (chrom, pos) = vcf1_data[0], vcf1_data[1]

        vcf1_ziptuple=zip(vcf1_samples, vcf1_data[9::] )
        vcf2_ziptuple=zip(vcf2_samples, vcf2_data[9::] )

        #filetr the ziptuple to contain only those samples common to both!
        filtered_vcf1= [x for x in vcf1_ziptuple if x[0] in common_samples]
        filtered__vcf2=  [x for x in vcf2_ziptuple if x[0] in common_samples]

        maf_vcf2 = calMaf(filtered__vcf2) #get the MAF for the site in the second VCF
        maf_vcf2= round ( maf_vcf2, 3)
        #print len(filtered_vcf1), len(filtered__vcf2)

        comparison_results = compare_genotypes(filtered_vcf1, filtered__vcf2)
      
        
        vcf1_formatfields=vcf1_formatstr.split(':')
        dp_index=vcf1_formatfields.index('DP')
        ra_index=vcf1_formatfields.index('RA')
        aa_index=vcf1_formatfields.index('AA')
        sr_index=vcf1_formatfields.index('SR')
        sa_index=vcf1_formatfields.index('SA')

        zero_counts_for_alt_homozref=0 #number of genotypes for which called homoz_ref by vcf_1 and het by vcf_2, no alt allele observed for vcf_1 genotype
        zero_counts_for_ref_homozalt=0 ##number of genotypes for which called homoz_alt by vcf_1 and het by vcf_2, no ref allele observed for vcf_2 genotype
        for i in range(0, len(comparison_results) ):
            (g1,g2,sample)=comparison_results[i]

            site_discordance[g1,g2]+=1

            if g2 == 1:
                if g1 == 1 or g1 == 1:
                #collect read depth info on site (ref and alt)
                    dp=int(getFormatfield(filtered_vcf1[i][1], dp_index) )
                    ra=int(getFormatfield(filtered_vcf1[i][1], ra_index) )
                    aa=int (getFormatfield(filtered_vcf1[i][1], aa_index) )
                    sr=( getFormatfield(filtered_vcf1[i][1], sr_index) )
                    sa= ( getFormatfield(filtered_vcf1[i][1], sa_index) )
                    
                #print g1, g2, filtered_vcf1[i], vcf1_formatstr, dp, ra, aa, sr, sa

                if g1 == 1:
                    if aa == 0 or aa==1:
                        zero_counts_for_alt_homozref+=1
                if g1 == 1:
                    if ra ==0 or ra==1:
                        zero_counts_for_ref_homozalt+=1
            

        #print site_discordance
        if site_discordance[1,1] !=0:
            nosample_refrate = round ( float(zero_counts_for_ref_homozalt)/float( site_discordance[1,1]),3 )
        else:
            nosample_refrate='NA'
        if site_discordance[1,1] !=0:
            nosample_altrate = round ( float(zero_counts_for_alt_homozref)/float( site_discordance[1,1]),3 )
        else:
            nosample_altrate='NA'
        
        vcf2_totalhetcalls= site_discordance[2,1] +  site_discordance[1,1] +  site_discordance[0,1]
        vcf1_hetdiscrepancy= float ( site_discordance[2,1] + site_discordance[0,1]) / float(vcf2_totalhetcalls)
        vcf1_hetdiscrepancy = round (vcf1_hetdiscrepancy, 3)

        vcf1_nocalls = site_discordance[3,0] + site_discordance[3,1] + site_discordance[3,2] + site_discordance[3,3]
        vcf2_nocalls = site_discordance[0,3] + site_discordance[1,3] + site_discordance[2,3] + site_discordance[3,3]


#        totalCalledGenotypes = computeCalledGenotypes(site_discordance)
#        totalChroms=2*totalCalledGenotypes

#        vcf1_altcount=computeRowMarginalAltCount(site_discordance)
#        vcf2_altcount =computeColMarginalAltCount(site_discordance)

        print chrom, pos, str(nosample_altrate), str( site_discordance[1,1] ), str(nosample_refrate),  str( site_discordance[1,1] ), str(vcf2_totalhetcalls), str(vcf1_hetdiscrepancy), maf_vcf2
        

if __name__ == "__main__":
    main()
