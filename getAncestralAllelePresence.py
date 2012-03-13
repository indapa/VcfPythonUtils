#!/usr/bin/env python
import sys
import os
import string
import re
from optparse import OptionParser
from VcfFile import *
""" get presence/absence of ancestral allele (AA) from VCF file """

def main():
    usage = "usage: %prog [options] arg"
    parser = OptionParser(usage)
    parser.add_option("--info", type="string", dest="infotag", help="INFO tag id that annotates what type of variant the VCF record is", default="TYPE")
    parser.add_option("--type", type="string", dest="variantype", help="type of variant (SNP INS DEL)", default='SNP')
    parser.add_option("--filter", type="string", dest="filter", help="extract records matching filter (default is PASS)", default='PASS')

    (options, args)=parser.parse_args()
    if options.infotag == "":
        sys.stderr.write("provide a value for --info parameter!\n")
        exit(1)
    if options.variantype == "":
        sys.stderr.write("provide a value of --type parameter!\n")
        exit(1)

    (options, args)=parser.parse_args()

    vcfilename=args[0]
    vcfh=open(vcfilename,'r')

    #instantiate a VcfFile object
    vcfobj=VcfFile(vcfilename)
    
    
    #parse its metainfo lines (ones that begin with ##)
    vcfobj.parseMetaLines(vcfh)
    vcfh.seek(0)


    #vcfobj.printMetaLines()
    vcfobj.parseHeaderLine(vcfh)
    #vcfobj.printHeaderLine()
    samplelist=vcfobj.getSampleList()
    newsamplestring='CHROM\tPOS\t'
    for s in samplelist:
        newsamplestring+= s+".1"+"\t"+s+".2\t"
    print newsamplestring
    descriptors = vcfobj.getMetaInfoDescription()
    
    infoids=[]
    for (tag, description) in descriptors:
        infoids.append(tag)

    if options.infotag not in infoids and options.infotag != 'QUAL':
        sys.stderr.write(options.infotag + " tag not in ##INFO headers!\n")
        exit(1)

    vcfh.seek(0)
    
    if options.variantype != None:
        pattern=options.infotag+'=('+options.variantype+')'

    aa_pattern='\tAA=([ACGTacgt])\t'

    ref_is_ancestral=0
    aa_pattern='(AA=([ACGTacgt]))'
    
  
    for dataline in vcfobj.yieldVcfDataLine(vcfh):
        fields=dataline.strip().split('\t')
        (chrom,pos,id,ref,alt,qual,filtercode,info)=fields[0:8]
            #check if its PASS filter
        if filtercode != options.filter and options.filter != None : continue

        ref_is_ancestral=0

        #check if its a SNP
        if re.search(pattern, info ) == None:
            continue
        
        if re.search(aa_pattern, info )== None:
            continue

        #if it gets here, record has an AA snp
        aa_allele=re.search(aa_pattern,info ).groups()[1]
        aa_allele=aa_allele.upper()
        #print chrom, pos, ref,aa_allele
        if aa_allele == ref:
            ref_is_ancestral = 1
        genotypes=[]
        genostrings=fields[9::]
        for gstring in genostrings:
            (gt, gq, gl)=gstring.split(':')
            genotypes.append(gt)
        #skip if not sites are called
        if '.' in genotypes:
            continue
        gzip=[]
        for g in genotypes:
            (a1,a2)=g.split('|')
            if ref_is_ancestral == 1:
               gzip.append( (a1,a2))
            else:
                gzip.append( ('0','0') )
        #print zip(samplelist, gzip)
        outstring=[chrom, pos]
        for ( sample, alleletuple ) in zip(samplelist, gzip):
            allestr="\t".join( list(alleletuple ))
            outstring.append(allestr)
        print "\t".join(outstring)



if __name__ == "__main__":
    main()
