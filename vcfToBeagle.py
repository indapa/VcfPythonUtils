#!/usr/bin/python
import sys
import os
import string
import re
from optparse import OptionParser
from vcfIO import *


def returnIdentifierLine(fh):
    """return Beagle id column headers with sample names taken from VCF filehandle """
    id_string=["I", "id"]
    samples=get_vcfsamples(fh)
    samples_diploid=[]
    
    for s in samples:
        samples_diploid.append(s)
        samples_diploid.append(s)

    samplestring= "\t".join(samples_diploid)
    idstring="I\tid\t"+samplestring
    return idstring



#beagle format described in http://faculty.washington.edu/browning/beagle/beagle_3.3_26Dec10.pdf
def main():
    """ create a Beagle v3.3 input file from a VCF file """
    usage = "usage: %prog [options] vcfile"
    parser = OptionParser(usage)

    (options, args)=parser.parse_args()

    vcfile=args[0]
    bglfile=vcfile.replace('.vcf', '.bgl')
    bgfh=open(bglfile,'w')

    try:
        sys.stderr.write("opening vcfile ....\n")
        vcf_fh=open(vcfile, 'r')
    except:
        sys.stderr.write("unable to open vcfile!\n")
            
    bglIdline=returnIdentifierLine(vcf_fh)
    bgfh.write(bglIdline+"\n")

    vcf_fh.seek(0)

    sys.stderr.write("writing Beagle file input file ...\n")
    #each t(uple) represents a snp locus with 2*N alleles
    for (chrom, pos, ref, alt,genotype_tuple) in get_vcftuples(vcf_fh):
        #print genotype_tuple
        alleles=[]
        for (sample, gt) in genotype_tuple:
            (allele1, allele2) = ('.', '.')
            strippedGenotype=stripGT(gt)
            #print strippedGenotype
            (allele1, allele2) = returnAlleles_unphased( strippedGenotype )
            
            if allele1 == '0':
                allele1=ref
            elif allele1 == '1':
                allele1 = alt
            elif allele1 == '.':
                allele1 = '?'
            else:
                pass

            if allele2 == '0':
                allele2=ref
            elif allele2 == '1':
                allele2=alt
            elif allele2 == '.':
                allele2='?'
            else:
                pass
            #print strippedGenotype, allele1, allele2, "\t",  ref, alt
            alleles.append(allele1)
            alleles.append(allele2)
        
        marker_id=chrom+"."+pos
        allele_str="\t".join(alleles)
        marker_line="\t".join(['M',marker_id, allele_str])
        bgfh.write(marker_line+"\n")

    bgfh.close()

if __name__ == "__main__":
    main()
