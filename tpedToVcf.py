#!/usr/bin/env python
import sys
import os
import string
import re
import bx.seq.twobit
from optparse import OptionParser
from common import grouper
from common import numericalGenotypes
import datetime


def printvcfHeader(tpedfile, tbfile):
    today=datetime.datetime.today()
    datestr=today.strftime("%m-%d-%y")
    #vcf_metainfolines
    vcf_metalines=[]
    vcf_metalines.append ( "##fileformat=VCFv4.1")
    vcf_metalines.append( "##fileDate="+datestr )
    vcf_metalines.append("##reference="+tbfile)
    vcf_metalines.append("##pedfile="+tpedfile)
    vcf_metalines.append("##INFO=<ID=NS,Number=1,Type=Integer,Description=\"Number of Samples With Data\">")
    vcf_metalines.append( "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">" )

    print "\n".join(vcf_metalines)
    

""" Given a tped/tfam file based on --file basename make a barebones VCF file
    A tped/tfam file is generated in PLINK, for more see here: http://pngu.mgh.harvard.edu/~purcell/plink/dataman.shtml#recode
    the only field in the format is GT
    No qual or filter value is given, only a '.'
"""
def main():
    usage = "usage: %prog [options] filebasename"
    parser = OptionParser(usage)
    parser.add_option("--file", type="string", dest="basename", help="basename of tped/tfam file")
    parser.add_option("--twobitfile", type="string", dest="tbf", help="2bit file of reference genome")


    (options, args)=parser.parse_args()
    

    try:
        sys.stderr.write("opening twobitfile...\n")
        twobit=bx.seq.twobit.TwoBitFile( open( options.tbf ) )
    except:
       sys.stderr.write("unable to open twobit file!\n")

    tfamfile=options.basename+".tfam"
    tpedfile=options.basename+".tped"

    tfamfh=open(tfamfile, 'r')
    samplenames=[]
    for line in tfamfh:
        (fid,iid,pid,mid,sex,pheno)=line.strip().split(' ')
        samplenames.append(iid)

    samplestring="\t".join(samplenames)


    tpedfh=open(tpedfile,'r')
    printvcfHeader(options.tbf, tpedfile)
    print "\t".join(["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER",  "INFO ",  "FORMAT", samplestring])
    for line in  tpedfh:
        fields=line.strip().split(' ')
        (chrom, snpid,cM,pos)=fields[0:4]
        start=int(pos)-1
        end=int(pos)
        
        try:
            sequence=twobit[chrom][start:end]
            sequence=sequence.upper()
        except:
            error="unable to fetch sequence from 2bit file!: " + chrom + " " + pos
            sys.stderr.write(error + "\n")
            exit(1)
        refbase=sequence
        #print chrom, pos,refbase
        genotypes=fields[4::]
        if len(genotypes)/2 != len(samplenames):
            sys.stderr.write("unequal numbers of genotypes and sample names!\n")
            sys.exit(1)
        observed_alleles=set(genotypes)
        altbases= list( observed_alleles - set(refbase) )
        alt='.'
        if len(altbases) == 0:
            alt='.'
        elif len(altbases) > 1:
            alt=",".join(altbases )
        else:
            alt=altbases[0]

        
        metainfo="\t".join([chrom,pos,snpid,refbase,alt,'.','.', 'NS='+str(len(samplenames)),'GT'])
        
        ngenotypes=[]
        for genotype in grouper(2, genotypes,'x'):
            genostr="".join(list(genotype) )
            ngenotypes.append( numericalGenotypes(refbase,alt, genostr) )

        #print genotypes
        #print ngenotypes
        goutput="\t".join(ngenotypes)
        print metainfo +"\t" +goutput



if __name__ == "__main__":
    main()
