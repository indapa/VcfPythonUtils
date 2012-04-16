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


        #if filtercode != filtercodeoption and filtercodeoption != None:
        #    continue


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
    parser.add_option("--filter", type="string", dest="filter", default=None, help="intersect records only set with filter (default is None")
    parser.add_option("--info", type="string", dest="infotag", help="INFO tag id that annotates what type of variant the VCF record is", default="TYPE")
    parser.add_option("--type", type="string", dest="variantype", help="type of variant (SNP INS DEL)", default="")
    parser.add_option("--noheader", action="store_true", dest="noheader", help="VCF file one  has no header line", default=False)

    (options, args)=parser.parse_args()

    sys.stderr.write("intersecting two files ...\n")
    
    vcf_file_one=args[0]
    in2_fname=args[1]

    

    if ".bed" in in2_fname:
        bitsets = binned_bitsets_from_file( open( in2_fname ) )

    if ".vcf" in in2_fname:
         bitsets = binned_bitsets_from_vcffile( in2_fname , options.filter)

   
    vcfobj=VcfFile(vcf_file_one)
    vcfh=open(vcf_file_one,'r')

    if options.noheader == False:
        vcfobj.parseMetaLines(vcfh)
        vcfobj.printMetaLines()
    




    descriptors = vcfobj.getMetaInfoDescription()
    infoids=[]
    for (tag, description) in descriptors:
        infoids.append(tag)

    if options.infotag  not in infoids and options.infotag != 'QUAL' and  options.infotag != "" and options.noheader == False:
        sys.stderr.write(options.infotag + " tag not in ##INFO headers!\n")
        exit(1)


    #vcfh.seek(0)

    #parse the header  line #CHROM and print it
    if options.noheader==False:
        vcfobj.parseHeaderLine(vcfh)
        vcfobj.printHeaderLine()
    
    if options.variantype != "":
        pattern=options.infotag+'=('+options.variantype+')'


    for dataline in vcfobj.yieldVcfDataLine(vcfh):
        fields=dataline.strip().split('\t')
        (chrom,pos,id,ref,alt,qual,filtercode,info)=fields[0:8]
        (start,end) = (int(pos)-1, int(pos))

        #pass the filter code
        if filtercode != options.filter and options.filter != None:
            continue

        #check to see if record is the correct variant TYPE
        if options.variantype != "":
            pattern=options.infotag+'=('+options.variantype+')'
            if re.search(pattern, info ) == None:
                continue
            else:
                value=re.search(pattern, info ).groups()[0]
                pass

        #chrom="chr"+chrom
        if chrom in bitsets and bitsets[chrom].count_range( start, end-start ) >= options.mincols:
            if not options.reverse:
                print dataline
        else:
            if options.reverse == True:
                print dataline
        


if __name__ == "__main__":
    main()
