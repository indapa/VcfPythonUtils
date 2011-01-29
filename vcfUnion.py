#!/usr/bin/python

from optparse import OptionParser
import sys
from bx.bitset import *
from bx.bitset_utils import *
from itertools import *
__author__="amit"
__date__ ="$Aug 19, 2010 7:15:06 PM$"


def splitConsecutive(chr,start,end):
    """
    bitset operation merge consecutive bits into a single bed interval
    but for our purposes we will split consecutive intervals into atomic 1 bp intervals
    """

    start = int(start)
    end = int(end)
    for i in range(start,end):
        print "%s\t%d\t%d" % ( chr, i, i+1 )

def binned_bitsets_from_vcffile( f, chrom_col=0, start_col=1,  upstream_pad=0, downstream_pad=0, lens={} ):
    """
    Read a file into a dictionary of bitsets. The defaults arguments

    - 'f' should be a file like object (or any iterable containing strings)
    - 'chrom_col', 'start_col', and 'end_col' must exist in each line.

    - if 'lens' is provided bitset sizes will be looked up from it, otherwise
      chromosomes will be assumed to be the maximum size
    """
    last_chrom = None
    last_bitset = None
    bitsets = dict()
    for line in f:
        if line.startswith("#") or line.isspace():
            continue
        fields = line.split('\t')
        
        strand = "+"
        chrom = fields[chrom_col]
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
    """ take the union of sites in an abritary number of VCF files and print the bed3 coordinates to STDOUT """
    """ uses the bx-python bitset API https://bitbucket.org/james_taylor/bx-python/wiki/Home """
    usage = "usage: %prog file1.vcf file2.vcf ....\nTakes the union of sites in an arbritary number of VCF files and prints the bed3 coordinates to STDOUT\n"
    parser = OptionParser(usage)

    (options, args) = parser.parse_args()
    vcf_files  = args[0:]
    #print vcf_files

    if vcf_files:
        input = chain( * imap( open, vcf_files ) )
    else:
        input = sys.stdin

    #print infile
    #vcfh = open(infile, 'r')
    bitsets = binned_bitsets_from_vcffile(input)
    

    for chrom in bitsets:
        bits = bitsets[chrom]
        end = 0
        while 1:
            start = bits.next_set( end )
            if start == bits.size: break
            end = bits.next_clear( start )
            if (end-start) >1:
                splitConsecutive(chrom, start, end)
            else:
                print "chr%s\t%d\t%d" % ( chrom, start, end )

   

if __name__ == "__main__":
    main()
