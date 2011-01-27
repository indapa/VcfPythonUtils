#!/usr/bin/python
import sys
import os
import string
import re
from optparse import OptionParser
from vcfIO import *


def printIdentifierLine(fh):
"""print Beagle id column headers with sample names taken from VCF filehandle """
    id_string=["I", "id"]
    samples=get_vcfsamples(fh)
    samples_diploid=[]
    
    for s in samples:
        samples_diploid.append(s)
        samples_diploid.append(s)

    samplestring= "\t".join(samples_diploid)

    print "I\tid\t",samplestring

#beagle format described in http://faculty.washington.edu/browning/beagle/beagle_3.3_26Dec10.pdf
def main():
""" create a Beagle v3.3 input file from a VCF file """
    usage = "usage: %prog [options] vcfile"
    parser = OptionParser(usage)

    (options, args)=parser.parse_args()

    vcfile=args[0]
    

    try:
        sys.stderr.write("opening vcfile ....\n")
        vcf_fh=open(vcfile, 'r')
    except:
        sys.stderr.write("unable to open vcfile!\n")
            
    printIdentifierLine(vcf_fh)
    vcf_fh.seek(0)
    
    #each t(uple) represents a snp locus with 2*N alleles
    for (chrom, pos, ref, alt,genotype_tuple) in get_vcftuples(vcf_fh):
        alleles=[]
        for (sample, gt) in genotype_tuple:
            (allele1, allele2) = ('.', '.')
            if stripGT(gt) == '.': # no genotype call represent as ? in beagle
                allele1='?'
                allele2='?'
            else:
               (allele1, allele2) = returnAlleles_unphased( stripGT(gt) )

            alleles.append(allele1)
            alleles.append(allele2)
        marker_id=chrom+"."+pos
        allele_str="\t".join(alleles)
        marker_line="\t".join(['M',marker_id, allele_str])
        print marker_line


if __name__ == "__main__":
    main()
