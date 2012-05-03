#!/usr/bin/env python
import sys
import os
import string
import re
from optparse import OptionParser

from VcfFile import *

def main():
    usage = "usage: %prog [options] file1.vcf file2.vcf"
    parser = OptionParser(usage)
    (options, args)=parser.parse_args()

    vcfileone=args[0]
    vcfiletwo=args[1]

    vcfh1=open(vcfileone, 'r')
    vcfh2=open(vcfiletwo, 'r')

    vcfobj1=VcfFile(vcfileone)
    vcfobj2=VcfFile(vcfiletwo)
    
    #parse its metainfo lines (ones that begin with ##)
    vcfobj1.parseMetaAndHeaderLines(vcfh1)
    vcfobj2.parseMetaAndHeaderLines(vcfh2)


    samples_fileone=vcfobj1.getSampleList()
    samples_filetwo=vcfobj1.getSampleList()
    print samples_fileone, samples_filetwo
    common_samples=set(samples_fileone).intersection( set(samples_filetwo) )
    N=len(common_samples)
    if N == 0:
        sys.stderr.write("no common samples between vcf files provided!")
        exit(1)



    """ kind of a convluluted way to go about iterating thru vcf datalines
        since the yield methods return a generator, get the geneator first
        then proceed to a while loop within a try/catch block, as the generators return
        an excption StopIteration when they are through"""

    vcf_generator1=vcfobj1.yieldVcfRecordwithGenotypes(vcfh1)
    vcf_generator2=vcfobj1.yieldVcfRecordwithGenotypes(vcfh2)
    try:
        while (1):
            vrec1= vcf_generator1.next()
            vrec2= vcf_generator2.next()

            print vrec1.toStringwithGenotypes()
            print vrec2.toStringwithGenotypes()

            zipListOne=vrec1.zipGenotypes(samples_fileone)
            zipListTwo=vrec1.zipGenotypes(samples_filetwo)

            zipListOne= [(samp, gt) for (samp, gt) in zipListOne if samp in common_samples ]
            zipListTwo= [ (samp, gt) for (samp, gt) in zipListTwo if samp in common_samples ]

            for i in range(0,N):
                print zipListOne[i][1].getAlleles()
                print zipListTwo[i][1].getAlleles()

                



    except StopIteration:
        sys.stderr.write("\ndone iterating through vcf records in both files\n")
    
    



if __name__ == "__main__":
    main()
