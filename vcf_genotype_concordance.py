#!/usr/bin/env python
import sys
import os
import string
import re
from optparse import OptionParser

import itertools
from VcfFile import *
import numpy as np


def typeofGenotype(allele1, allele2):
    """ I really should be a python version of a typedef here, but dont know how
        hom_ref =0 het =1 hom_nonref=2 no_call=3                              """

    if allele1 == '0' and allele2 == '0': return 0

    if allele1 == '0' and allele2== '1': return 1
    if allele1 =='1' and allele2 == '0': return 1

    if allele1== '1' and allele2== '1': return 2

    if allele1== '.' or allele2 == '.': return 3
    if allele1 == '.' and allele2 == '.': return 3

def calMaf (g):
    """ given a list in the form [ (sample, VcfGenotype obj) ] determine the minor allele freq= alt/2*called genotypes """
    total=0
    alt_count=0
    for (sample, genotype) in g:
        (p1,p2) = genotype.getAlleles()
        
        alleletype=typeofGenotype(p1,p2)

        if alleletype == 3:
            continue
        total+=1
        if alleletype ==1 : 
            alt_count+=1
        elif alleletype == 2:
            alt_count +=2
        else:
            pass

    if float(2*total) == 0:
        return (round(-1,3))
    
    maf = float(alt_count) / float(2*total)
    #print float(2*total)
    #print p1, p2, sample, maf, alt_count, total
    return  ( round(maf,3) )

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
    #print gtm

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
        nrd_homoz_nonref = 'NA'
    else:
        nrd_homoz_nonref= round ( float(discord_homoz_nonref)/float(total_homoz_nonref), 3 )


    return (nrd_homoz_ref, nrd_het, nrd_homoz_nonref,
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



def compare_genotypes(g1, g2):
    """ given two lists of the form [ (sample, VcfGenotype obj), ....  ], iterate thru and compare genotypes per individual    """
    """ return list of [ ( g1_alleletype, g2_alleletype, samplename) ] where alleletype is [0, homref; 1, het; 2 hom_nonref; 3, nocall """
    compared_genotypes=[] # [ (g1_type, g2_type, sample) ... ]

    if len(g1) != len(g2):
        sys.stderr.write("cannot compare phased unphased genotypes; unequal size of lists!")
    else:
        for i in range( 0, len( g1 ) ):
            if g1[i][0] != g2[i][0]:
                sys.stderr.write("samples don't match in genotypes comparison!")
                exit(1)
            (p1,p2) = g1[i][1].getAlleles()
            (u1, u2) = g2[i][1].getAlleles()

            g1_alleletype=typeofGenotype(p1,p2)
            g2_alleletype=typeofGenotype(u1,u2)
            
            compared_genotypes.append( (g1_alleletype, g2_alleletype, g1[i][0])  )

    return  compared_genotypes

"""   generate concordance metrics accocrding to this given two vcf files with the same data, but assigned/called variants
      from different methods
      http://www.broadinstitute.org/gsa/wiki/index.php/GenotypeConcordance """

def main():
    usage = "usage: %prog [options] file1.vcf file2.vcf"
    parser = OptionParser(usage)
    parser.add_option("--filter", type="string", dest="filter", default=None, help="analyze records only set with filter (default is None")
    parser.add_option("--info", type="string", dest="infotag", help="INFO tag id that annotates what type of variant the VCF record is", default="TYPE")
    parser.add_option("--type", type="string", dest="variantype", help="type of variant (SNP INS DEL)", default=None)

    (options, args)=parser.parse_args()

    sitefh = open('site.nrd.nrs.txt', 'w')
   
    vcfileone=args[0]
    vcfiletwo=args[1]

    headerline="\t".join(['#chrom', 'pos', 'ref','alt', 'total_samples_vcf1', 'total_called_genotypes_vcf1',
                               '#chromosomes','AC_'+vcfileone,'noCalls_'+vcfileone,'AC_'+vcfiletwo,'noCalls_'+vcfiletwo,'maf_'+vcfiletwo,
                                'site_NRS','site_NRD', 'site_nrd_homoz_ref','site_nrd_het', 'site_nrd_homoz_nonref'])
    sitefh.write(headerline+"\n")


    vcfh1=open(vcfileone, 'r')
    vcfh2=open(vcfiletwo, 'r')

    vcfobj1=VcfFile(vcfileone)
    vcfobj2=VcfFile(vcfiletwo)
    
    #parse its metainfo lines (ones that begin with ##)
    vcfobj1.parseMetaAndHeaderLines(vcfh1)
    vcfobj2.parseMetaAndHeaderLines(vcfh2)


    samples_fileone=vcfobj1.getSampleList()
    samples_filetwo=vcfobj2.getSampleList()
    sys.stderr.write("common samples between  files:\n" + "file one \t".join(samples_fileone) + "\nfile two\n" + "\t".join(samples_filetwo) + "\n")
    common_samples=set(samples_fileone).intersection( set(samples_filetwo) )
    N=len(common_samples)

    if N == 0:
        sys.stderr.write("no common samples between vcf files provided!")
        exit(1)

    sitetotalAA=0
    sitetotalAA_discord=0

    sitetotalAB=0
    sitetotalAB_discord=0

    sitetotalBB=0
    sitetotalBB_discord=0



    """ a discordance table per sample """
    discordance_dict={}
    for s in common_samples:
        discordance_dict[s]=np.matrix( [ [ 0,0,0,0 ], [ 0,0,0,0 ], [ 0,0,0,0 ], [ 0,0,0,0 ] ] )


    """ Make an iterator that aggregates elements from each of the iterables.
    Like zip() except that it returns an iterator instead of a list.
    Used for lock-step iteration over several iterables at a time
    http://docs.python.org/library/itertools.html#itertools.izip"""

    vcf_gen1=vcfobj1.yieldVcfRecordwithGenotypes(vcfh1)
    vcf_gen2=vcfobj1.yieldVcfRecordwithGenotypes(vcfh2)

    
    sys.stderr.write("computing non-ref discrepancy and non-refernce sensititvity\n")
    for vrec1, vrec2 in itertools.izip(vcf_gen1, vcf_gen2):
            #print vrec1.toStringwithGenotypes()
            #print vrec2.toStringwithGenotypes()

            if vrec1.getRef() != vrec2.getRef():
                sys.stderr.write("Vcf records don't match in chrom!\n")
                outstring="\t".join([ vrec1.getChrom(), vrec1.getPos(), vrec2.getChrom(), vrec2.getPos() ])
                sys.stderr.write(outstring+"\n")
                #continue
            if vrec1.getPos() != vrec2.getPos():
                sys.stderr.write("Vcf records don't match in positions!\n")
                outstring="\t".join([ vrec1.getChrom(), vrec1.getPos(), vrec2.getChrom(), vrec2.getPos() ])
                sys.stderr.write(outstring+"\n")
                #continue

            if vrec1.getRef() != vrec2.getRef():
                sys.stderr.write("Vcf records don't match in reference allele!\n")
                outstring="\t".join([ vrec1.getChrom(), vrec1.getPos(), vrec1.getRef(), vrec1.getAlt(),  vrec2.getChrom(), vrec2.getPos(), vrec2.getRef(), vrec2.getAlt()  ])
                sys.stderr.write(outstring+"\n")
                #continue
            if vrec1.getAlt() != vrec2.getAlt():
                sys.stderr.write("Vcf records don't match in alt allele!\n")
                outstring="\t".join([ vrec1.getChrom(), vrec1.getPos(), vrec1.getRef(), vrec1.getAlt(),  vrec2.getChrom(), vrec2.getPos(), vrec2.getRef(), vrec2.getAlt()  ])
                sys.stderr.write(outstring+"\n")
                #continue

            vrec1info=vrec1.getInfo()
            vrec2info=vrec2.getInfo()

            """ check to see vrecs are the correct variant type
                and to see that vrec1 is of same variant type as vrec2
                so we compare apples to apples                    """
            if options.variantype != None:
                pattern=options.infotag+'=('+options.variantype+')'
                if re.search(pattern, vrec1info ) == None or re.search(pattern, vrec2info) == None:
                    continue
                vrec1type=re.search(pattern, vrec1info ).groups()[0]
                vrec2type=re.search(pattern, vrec2info ).groups()[0]
                if vrec1type != vrec2type:
                   sys.stderr.write("Vcf records don't match in variant type!\n")
                   outstring="\t".join([ vrec1.getChrom(), vrec1.getPos(),vrec1type, vrec2.getChrom(), vrec2.getPos(), vrec2type ])
                   sys.stderr.write(outstring+"\n")
                #print vrec1type, vrec2type
            else:
                pattern=options.infotag+'=(\w+)'
                vrec1type=re.search(pattern, vrec1info ).groups()[0]
                vrec2type=re.search(pattern, vrec2info ).groups()[0]
                if vrec1type != vrec2type:
                    sys.stderr.write("Vcf records don't match in variant type!\n")
                    outstring="\t".join([ vrec1.getChrom(), vrec1.getPos(),vrec1type, vrec2.getChrom(), vrec2.getPos(), vrec2type ])
                    sys.stderr.write(outstring+"\n")
                #print vrec1type, vrec2type

            """ list of tuples [ (sample, genotype object), .... ] """
            vrec1_ziptuple=vrec1.zipGenotypes(samples_fileone)
            vrec2_ziptuple=vrec2.zipGenotypes(samples_filetwo)
            
            """filter the list of tuples to contain only those samples common to both!"""
            filtered_vrec1= [x for x in vrec1_ziptuple if x[0] in common_samples]
            filtered__vrec2=  [x for x in vrec2_ziptuple if x[0] in common_samples]

            maf_vrec2 = calMaf(filtered__vrec2) #get the MAF for the site in the second VCF

            """the genotype comparison matrix for the VCF site
            reset everytime a new site is analyzed"""
            site_discordance=np.matrix( [ [ 0,0,0,0 ], [ 0,0,0,0 ], [ 0,0,0,0 ], [ 0,0,0,0 ] ] )
            #print site_discordance
            """ the following returns a list of tuples in the following way:
                [(vrec1_alleletype,vrec2_alleletype, sample)]
                allele_type is a number between 0-3 where 0:homozref 1:het 2:homozref"""
            comparison_results=compare_genotypes(filtered_vrec1, filtered__vrec2)
            #print comparison_results

            #this is where we increment the counts in the 
            # GenoTypeMatrix (gtm) gets updated here for the site and sample level
            for (g1, g2, sample) in comparison_results: #iterate thru the compariosn results
                discordance_dict[sample][g1,g2]+=1      #incrment the per-sample counts of genotype compariosns results
                site_discordance[g1,g2]+=1             #incrment the per-site counts of genotype compariosns results
            
            site_nrs=computeNRS( site_discordance )
            site_nrd=computeNRD ( site_discordance )
            (site_nrd_homoz_ref, site_nrd_het, site_nrd_homoz_nonref, aa, ab, bb )=computeNRD_class( site_discordance)
            sitetotalAA+=aa[0]
            sitetotalAA_discord+=aa[1]

            sitetotalAB+=ab[0]
            sitetotalAB_discord+=ab[1]

            sitetotalBB+=bb[0]
            sitetotalBB_discord+=bb[1]

            #print site_discordance
            #print site_nrs
            #print site_nrd

            vrec1_nocalls = site_discordance[3,0] + site_discordance[3,1] + site_discordance[3,2] + site_discordance[3,3]
            vrec2_nocalls = site_discordance[0,3] + site_discordance[1,3] + site_discordance[2,3] + site_discordance[3,3]


            totalCalledGenotypes = computeCalledGenotypes(site_discordance)
            totalChroms=2*totalCalledGenotypes

            vrec1_altcount=computeRowMarginalAltCount(site_discordance)
            vrec2_altcount =computeColMarginalAltCount(site_discordance)

            outstr= "\t".join( [ vrec1.getChrom(), vrec1.getPos(), vrec1.getRef(),
                             vrec1.getAlt() , str( len(filtered_vrec1) ),
                             str(totalCalledGenotypes), str(totalChroms), str(vrec1_altcount),
                             str(vrec1_nocalls), str(vrec2_altcount), str(vrec2_nocalls) ])

            sitefh.write(outstr +"\t")
            sitefh.write('%.2f'%maf_vrec2+'\t')
            sitefh.write( str(site_nrs)+'\t')
            sitefh.write( str(site_nrd)+'\t')
            sitefh.write( str(site_nrd_homoz_ref) + '\t')
            sitefh.write( str(site_nrd_het) + '\t')
            sitefh.write( str(site_nrd_homoz_nonref) + '\n')

            #break

    sitefh.close()
    #print "---"
    """broke out of the for vrec1, vrec2 in itertools.izip(vcf_gen1, vcf_gen2) loop
     so we finished comparing the two vcf files
     now collect and write per sample nrd and nrs results """

    outfh=open('sample.nrd.nrs.txt', 'w')
    headerstring = "\t".join( [ "sample","NRD", "NRS", "totalGenotypes", "noCallsEval", "noCallsComparison"] )
    outfh.write(headerstring+"\n")
    for sample in discordance_dict.keys():
        #print sample
        #print discordance_dict[sample]
        #print "total gneotypes: ", np.sum( discordance_dict[sample] )

        nrc=computeNRS( discordance_dict[sample] )
        nrd=computeNRD(  discordance_dict[sample]  )
        #print nrc
        #print nrd
        gtm= discordance_dict[sample]
        print sample
        print gtm
        print "\n"
        eval_nocalls = gtm[3,0] + gtm[3,1] + gtm[3,2] + gtm[3,3]
        comparison_nocalls = gtm[0,3] + gtm[1,3] + gtm[2,3] + gtm[3,3]
        outstring = "\t".join( [ sample, str(nrd), str(nrc),  str ( np.sum( discordance_dict[sample] ) ) , str(eval_nocalls) , str(comparison_nocalls ) ] )
        outfh.write(outstring+"\n")

    outfh.close()
if __name__ == "__main__":
    main()
