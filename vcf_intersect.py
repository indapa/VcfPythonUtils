#!/usr/bin/python
import sys
import os
import string
import re
from optparse import OptionParser
from VcfFile import *
from bx.bitset import *
from bx.bitset_builders import *


def binned_bitsets_from_vcffile( vcfilename, chrom_col=0, start_col=1,  upstream_pad=0, downstream_pad=0, lens={} ):
    """
    Read a vcffile into a dictionary of bitsets. The defaults arguments

    - 'vcfilename' should be a filename for vcf file
    - 'chrom_col', 'start_col', and 'end_col' must exist in each line.

    - if 'lens' is provided bitset sizes will be looked up from it, otherwise
      chromosomes will be assumed to be the maximum size

    - the bitset interval made into a   zero-based, half-open interval!!!!!!!

    """
    last_chrom = None
    last_bitset = None
    bitsets = dict()
    MAX=2147483647

    vcfobj=VcfFile(vcfilename)
    fh=open(vcfilename,'r')

    for vrec in vcfobj.yieldVcfRecord(fh):

        filtercode = vrec.getFilter()
        chrom = vrec.getChrom()
        pos=int( vrec.getPos() )
        
        if filtercode != 'PASS':
            if filtercode == '.':
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
        start, end = (pos-1, pos)

        if upstream_pad: start = max( 0, start - upstream_pad )
        if downstream_pad: end = min( size, end + downstream_pad )
        if start > end: warn( "Interval start after end!" )
        last_bitset.set_range( start, end-start )
    fh.close()
    return bitsets

def main():
    
    usage = "usage: %prog [options] vcf_file_one vcf|bed_file_two\n\nFind regions in the first vcf file that overlap regions of the second vcf or bed file\n"
    parser = OptionParser(usage)
    parser.add_option("--minCols", type="int", dest="mincols", default=1, help="mininum basepair overlap (default is one)")
    parser.add_option("--v", action="store_true", dest="reverse",  help="Print regions in first vcf  that DO NOT overlap second vcf|bed file")
    #parser.add_option("--sites", action="store_true", dest="sites", help="print only site information")
    (options, args)=parser.parse_args()

    sys.stderr.write("intersecting two files ...\n")
    
    vcf_file_one=args[0]
    in2_fname=args[1]

    

    if ".bed" in in2_fname:
        bitsets = binned_bitsets_from_file( open( in2_fname ) )

    if ".vcf" in in2_fname:
         bitsets = binned_bitsets_from_vcffile( in2_fname )

   
    vcfobj=VcfFile(vcf_file_one)
    vcfh=open(vcf_file_one,'r')
    vcfobj.parseMetaLines(vcfh)

    #print meta INFO, FORMAT, and FILTER lines
    vcfobj.printMetaLines()

    vcfh.seek(0)

    #parse the header  line #CHROM and print it
    vcfobj.parseHeaderLine(vcfh)
    vcfobj.printHeaderLine()
    
    for dataline in vcfobj.yieldVcfDataLine(vcfh):
        fields=dataline.strip().split('\t')
        (chrom,pos,id,ref,alt,qual,filtercode,info)=fields[0:8]
        (start,end) = (int(pos)-1, int(pos))
        if filtercode != 'PASS':
            continue
        chrom="chr"+chrom
        if chrom in bitsets and bitsets[chrom].count_range( start, end-start ) >= options.mincols:
            if not options.reverse:
                print dataline
        else:
            if options.reverse == True:
                print dataline
        


if __name__ == "__main__":
    main()
