#!/usr/bin/python
import sys
import os
import string
import re
from optparse import OptionParser
from vcfIO import *

def main():
    """ given a VCF file determine sample callrate and snp callrate. Results written to snpcallrate.txt and samplecallrate.txt """
    
    usage = "usage: %prog [options] vcf_file"
    parser = OptionParser(usage)


    (options, args)=parser.parse_args()

    vcfile=args[0]


    try:
        vcf_fh=open(vcfile, 'r')
    except:
        sys.stderr.write("unable to open vcfile!\n")

    try:
        snpcallrate_fh=open("snpcallrate.txt", 'w')
    except:
        sys.stderr.write("cannot create file snpcallrate.txt")
    
    try:
        samplecallrate_fh=open("samplecallrate.txt", 'w')
    except:
        sys.stderr.write("cannot create file samplecallrate.txt")

    samplecallrate_fh.write("sample\tcalled\ttotal\tcallrate\n")
    snpcallrate_fh.write("chr\tpos\tcalled\ttotal\tcallrate\n")
    
    samples=get_vcfsamples(vcf_fh)
    total_samples=len(samples)
    #print total_samples
    sampleCalls={} #key sample name value #called genotypes
    for s in samples: sampleCalls[s]=0
    
    vcf_fh.seek(0)

    totalSnps=0
    #snp callrate: total samples with called genotypes/total number of samples
    #sample callrate: total snp calls for indiv/total number of snps in VCF

    for (chrom, pos, ref, alt,genotype_tuple) in get_vcftuples(vcf_fh):
        totalSnps+=1
        calledsamples=0 #denominator for snp callrate
        for (sample, gt) in genotype_tuple:
            if "."  not in stripGT(gt):
                sampleCalls[sample]+=1
                calledsamples+=1
        outstring="\t".join([chrom,pos, str(calledsamples), str(total_samples), str(float(calledsamples)/float(total_samples)) ])
        snpcallrate_fh.write(outstring+"\n")
    #sample callrate
    for s in sampleCalls.keys():
        callrate=float(sampleCalls[s])/float(totalSnps)
        outstring="\t".join([s, str(sampleCalls[s]), str(totalSnps), str(callrate)])
        samplecallrate_fh.write(outstring+"\n")

if __name__ == "__main__":
    main()
