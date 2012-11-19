#!/usr/bin/python
import sys
import os
import string
import re
from optparse import OptionParser
from bx.align import maf
from vcfIO import *
import subprocess


def build_index(in_file):
    """ index a maf file """
    cl = ["maf_build_index.py", in_file]
    subprocess.check_call(cl)

def get_maf_slice ( maf_index,  org, chrom , start, end, assembly):
    
    region_name = "%s.%s" % (org, chrom)

    assembly_seq='.' # sequence of assembly you want to extract that is aligned against org

    for align in maf_index.get(region_name, start, end):
        region_align = align.slice_by_component(region_name, start, end)
        
        seqs_by_org = dict()
        for component in region_align.components:
            seqs_by_org[component.src] = component.text
        
        orgs=seqs_by_org.keys()

        for org in orgs:
            #print org
            (assembly_id, chromosome) = org.split('.',1)
            if assembly_id == assembly:
                assembly_seq=seqs_by_org[org]
    assembly_seq=assembly_seq.replace('-','')
    return assembly_seq.upper()

def main():
    """ extract  bases from MAF alignment  given positions out of a VCF file
        based on code found here: http://biostar.stackexchange.com/questions/3703/help-with-mafs-in-bx-python """
    usage = "\nusage: %prog [options] file.vcf\n Extract aligning base from  assembly (--assembly)  from MAF with aligning organism (--org) given a coordinate from VCF file. More info on MAF format here - http://genome.ucsc.edu/FAQ/FAQformat.html#format5   "
    parser = OptionParser(usage)
    parser.add_option("--maf", type="string", dest="maf", help="maf file")
    parser.add_option("--chrom", type="string", dest="chrom", help="chrom")
    parser.add_option("--org", type="string", dest="org", default = "hg19", help="organism (default hg19")
    parser.add_option("--assembly", type="string", dest="assembly", help="name of assembly whose aligning base (if any) against org you want to extract from MAF")
    (options, args)=parser.parse_args()


    vcf_file=args[0]
    try:
        vcf_fh=open(vcf_file, 'r')
    except:
        sys.stderr.write("cannot open vcf file!")

    
    #build maf index if it doesn't exist already
    index_file = options.maf + ".index"
    #print index_file
    if not os.path.exists(index_file):
        sys.stderr.write("indexing maf...\n")
        build_index(options.maf)
    maf_index = maf.Indexed(options.maf, index_file)


    headerline="\t".join( [ '#CHROM', 'POS',   'REF', 'ALT', options.assembly, 'FILTER'])
    print headerline
    sys.stderr.write("iterating thru VCF and parsing MAF....\n")
    for vcf_line in get_vcfdataline_passfilter(vcf_fh):
        vcf_data=split_vcfdataline(vcf_line)
        (chr, pos) = vcf_data[0:2]
        (id, ref, alt, qual,filter) = vcf_data[2:7]
        if filter=='.': filter="PASS"
        chr="chr" + chr
        if chr != options.chrom: continue
        start=int(pos)-1
        end = int(pos)
        
        assembly_seq= get_maf_slice ( maf_index, options.org, options.chrom, start, end, options.assembly)
        outstring= "\t".join( [ chr, pos, ref, alt,  assembly_seq, filter ] )
        print outstring

if __name__ == "__main__":
    main()
