#!/usr/bin/python


"""
using the bx-python interval API compare genotypes in the second VCF in the ones in the first VCF
Note, VCFs don't have to have the same number of individuals. For each site in the second VCF that overlaps the first,
only genotypes of samples common to both are compared.

Rather than use the bx-python bx.bitset  API for interval operation(s), I use the  bx.intervals because I can keep hold of annotations
of intervals.

"""


import sys
from optparse import OptionParser
from bx.intervals import *
from vcfIO import *

def makeVcfInterval(vcffh):

    sys.stderr.write("making vcf intervals ...\n")
    vcf_ranges = {}
    names=[]
    samplenames = get_vcfsamples(vcffh)

    for line in vcffh:
        fields = split_vcfdataline( line.strip() )
        #print fields


        genotypes = fields[9::]
        chrom,start,id,ref,alt,qual,filter,info = fields[0], int(fields[1]), fields[2], fields[3], fields[4], fields[5], fields[6], fields[7]
        #print filter
        if filter != 'PASS':
            if filter == '.':
                pass
            else:
                continue


        
        #otherstring = "\t".join(fields[9:])
        sample_genotype = zip( samplenames, fields[9::] )

        if chrom not in vcf_ranges: vcf_ranges[chrom] = Intersecter()

        vcf_ranges[chrom].add_interval( Interval( start-1, start,strand=1, name=id, info={'chr':chrom, 'ref':ref, 'alt':alt, 'qual':qual, 'filter':filter, 'info':info, 'samples':names, 'genotypes':genotypes , 'genotype':sample_genotype} ) )
        
    #print vcf_ranges
    return  vcf_ranges




def main():
    description="using the bx-python interval API compare genotypes in the second VCF in the onesin the first VCF.\nNote, VCFs don't have to have the same number of individuals.\n For each site in the second VCF that overlaps the first,only genotypes of samples common to both are compared. "
    usage = "\nusage: %prog [options] vcf_fileone\n" + "\n"+description
    parser = OptionParser(usage)
    #parser.add_option("--name", type="string|int|boolean", dest="name", help="help string")

    (options, args)=parser.parse_args()
    if len(args) !=2:
        sys.stderr.write("please provide 2 VCF file!\n")
        print usage
        exit(1)
    (in_fname1, in_fname2) = args[0:2]

    vcf_fh1=open(in_fname1, 'r')
    vcf_fh2=open(in_fname2, 'r')

    vcf_iname2_intervals=makeVcfInterval(vcf_fh2)

    vcf_fh1_samples = get_vcfsamples( vcf_fh1 ) #samples from the first vcf file
    #print vcf_fh1_samples
    
    for datafields in get_vcfdataline_passfilter( vcf_fh1 ):
        chr, start = datafields[0], int( datafields[1] )
        genotypes = datafields[9::]
        #ziptuple= zip( vcf_fh1_samples, genotypes )
        #print chr, start, ziptuple
        if chr not in vcf_iname2_intervals.keys():
            continue
        if len( vcf_iname2_intervals[chr].find( start-1, start ) ) > 0:
            print chr, start, datafields

if __name__ == "__main__":
    main()
