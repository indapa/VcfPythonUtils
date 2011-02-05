#!/usr/bin/python
"""
Find regions in the first vcf file that overlap regions of the second vcf or bed file
Similiar to bed_intersect.py from bx-python: https://bitbucket.org/james_taylor/bx-python/src/14b6a6c95da6/scripts/bed_intersect.py

Requires bx-python to be installed as it uses the bitset and bitet_builders api

VCF spec http://www.1000genomes.org/wiki/Analysis/Variant%20Call%20Format/vcf-variant-call-format-version-40
doesnt seem to require 'chr' prefix but most bedfiles I use have the 'chr' prefix in the target name, so I add the 'chr' when parsing the chrom field in a VCF

output is written to STDOUT

"""


from bx.bitset import *
from bx.bitset_builders import *
from vcfIO import *

from optparse import OptionParser

def binned_bitsets_from_vcffile( f, chrom_col=0, start_col=1,  upstream_pad=0, downstream_pad=0, lens={} ):
    """
    Read a vcffile into a dictionary of bitsets. The defaults arguments

    - 'f' should be a file like object (or any iterable containing strings)
    - 'chrom_col', 'start_col', and 'end_col' must exist in each line.

    - if 'lens' is provided bitset sizes will be looked up from it, otherwise
      chromosomes will be assumed to be the maximum size

    - the bitset interval made into a   zero-based, half-open interval!!!!!!!

    """
    last_chrom = None
    last_bitset = None
    bitsets = dict()
    MAX=2147483647

    for line in f:
        if line.startswith("#") or line.isspace():
            continue
        fields = line.split('\t')

        strand = "+"
        chrom = fields[chrom_col]
        filter  = fields[6]
       
        if filter != 'PASS':
            if filter == '.':
                pass
            else:
                continue


        chrom="chr"+chrom
        if chrom != last_chrom:
            if chrom not in bitsets:
                if chrom in lens:
                    size = lens[chrom]
                else:
                    size = MAX
                bitsets[chrom] = BinnedBitSet( size )
            last_chrom = chrom
            last_bitset = bitsets[chrom]
        start, end = int( fields[start_col])-1 , int( fields[start_col] )

        if upstream_pad: start = max( 0, start - upstream_pad )
        if downstream_pad: end = min( size, end + downstream_pad )
        if start > end: warn( "Interval start after end!" )
        last_bitset.set_range( start, end-start )
    return bitsets


def main():
    usage = "usage: %prog [options] vcf_file_one vcf|bed_file_two\n\nFind regions in the first vcf file that overlap regions of the second vcf or bed file\n"
    parser = OptionParser(usage)
    parser.add_option("--minCols", type="int", dest="mincols", default=1, help="mininum basepair overlap (default is one)")
    (options, args)=parser.parse_args()
    
    vcf_file_one=args[0]
    in2_fname=args[1]

    vcf_fh=open(vcf_file_one, 'r')

    if ".bed" in in2_fname:
        bitsets = binned_bitsets_from_file( open( in2_fname ) )

    if ".vcf" in in2_fname:
         bitsets = binned_bitsets_from_vcffile( open(in2_fname) )

    for dataline in vcf_fh:
        if '#' in dataline:
            print dataline.strip()
            continue
        fields=dataline.split('\t')
        (chrom,pos) =fields[0:2]
        start=int(pos)-1
        end=int(pos)
        #print chr, str(start), str(end)
        chrom="chr"+chrom
        if chrom in bitsets and bitsets[chrom].count_range( start, end-start ) >= options.mincols:
            print dataline.strip()




if __name__ == "__main__":
    main()
