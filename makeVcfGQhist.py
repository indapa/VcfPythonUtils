#!/usr/bin/python2.6
import sys
import os
import string
import re
from optparse import OptionParser
from vcfIO import *
import numpy as np
import matplotlib.pyplot as plt

""" given a VCF file with GQ in the format string plot a histogram of genotype qualties with help from NumPy and matplotlib """


def main():
    usage = "usage: %prog file.vcf\nGiven a VCF file with GQ in the format string\n plot a histogram of genotype qualties"
    parser = OptionParser(usage)
    

    (options, args)=parser.parse_args()

    vcfile = args[0]

    try:
        vcf_fh=open(vcfile, 'r')
    except:
        sys.stderr.write("cannot open vcf file!\n")
        exit(1)

    vcf_samples=get_vcfsamples(vcf_fh)
    vcf_fh.seek(0)

    genotype_quals=[]

    while 1:
        if '#CHROM' in vcf_fh.readline(): break

    sys.stderr.write("collecing genotype qualities ...\n")
    while(1):
        vcf_line = vcf_fh.readline()
        if vcf_line == '': break
        
        vcf_data=split_vcfdataline(vcf_line)
        formatstr=vcf_data[8]

        vcf_ziptuple=zip(vcf_samples, vcf_data[9::] )
        #for x in vcf_ziptuple:
        #    print x
        for gq in yieldGQ(vcf_ziptuple, formatstr):
            genotype_quals.append(int(gq) )


    sys.stderr.write("writing hisogram to GQhist.png...\n")
    fig = plt.figure()
    ax = fig.add_subplot(111)
    pdf, bins, patches = ax.hist(np.array ( genotype_quals ), 50, normed=1 )
   
    ax.set_xlabel('Genotype Quality')
    ax.set_ylabel('Probability')
    ax.grid(True)
    plt.title('histogram of genotype qualities')
    plt.savefig('GQhist.png')
    



if __name__ == "__main__":
    main()
