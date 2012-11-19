#!/usr/bin/python
import sys
import os
import string
import re
from optparse import OptionParser
from VcfFile import *
def main():
    usage = "usage: %prog [options] arg"
    parser = OptionParser(usage)
    

    (options, args)=parser.parse_args()

    vcfilename=args[0]
    
    vcfh=open(vcfilename, 'r')

    #instantiate a VcfFile object
    vcfobj=VcfFile(vcfilename)
    #parse its metainfo lines (ones that begin with ##)
    vcfobj.parseMetaLines(vcfh)
    #vcfobj.printMetaLines()
    vcfh.seek(0)
    vcfobj.parseHeaderLine(vcfh)
    samplelist=vcfobj.getSampleList()
    samplestrng="\t".join(samplelist)
    #vcfobj.printHeaderLine()
    headerstring="\t".join(['CHROM', 'POS', 'ID', 'REF', 'ALT', 'GENE', 'GENCODE_ID', 'STRAND', 'TYPE'])
    print headerstring, "\t", samplestrng
    for dataline in vcfobj.yieldVcfDataLine(vcfh):
        GT=[]
        #print dataline.strip()
        fields=dataline.strip().split('\t')
        (chrom,pos,id,ref,alt,qual,filtercode,info)=fields[0:8]
        genotypes=fields[9::]
        for g in genotypes:
            gfields=g.split(':')
            GT.append(gfields[0])
        
        infofields=info.split(';')
        VAT=infofields[-1]
        VATfields=VAT.split(':')
        #print VATfields
        (va, gene, gencode,strand, type)=VATfields[0:5]
        #print va, gene, gencode, strand, type
        gstring = "\t".join(GT)
        outstring = "\t".join( [ chrom, pos, id, ref, alt, gene, gencode, strand, type ] )
        print outstring , "\t", gstring



if __name__ == "__main__":
    main()
