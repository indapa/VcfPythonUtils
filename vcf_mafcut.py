#!/usr/bin/python
import sys
import os
import string
import re
from optparse import OptionParser

from VcfFile import *

""" extract out vcf data lines that meet minor allele freq thresholds """
""" for a MAF to be  extracted there needs to be a tag in the INFO field correpsonding to allele frequency """
def main():
    usage = "usage: %prog [options] maf file.vcf"
    parser = OptionParser(usage)
  
    parser.add_option("--maftag", type="string", dest="maftag", help="INFO tag id that annotates the allele freq of the record", default="AF")
    parser.add_option("--variantag", type="string", dest="vtag", help="INFO tag that annotates the type of variant type", default="VT")
    parser.add_option("--variantype", type="string", dest="variantype", help="type of variant (SNP INS DEL)", default=None)
    parser.add_option("--filter", type="string", dest="filter", help="extract records matching filter (default is None)", default=None)
    parser.add_option("--noheader", action="store_true", dest="noheader", help="VCF file  has no header file", default=False)
    parser.add_option("--quiet", action="store_true", dest="quiet", help="don't print vcf output to stdout", default=False)
    parser.add_option("--leq", type="float", dest="leq", default=1.0, help="keep variants with AF <= (default 1)")
    parser.add_option("--geq", type="float", dest="geq", default=0.0, help="keep variants with AF >= (default 0)")
    (options, args)=parser.parse_args()

    

    if len(args)!=1:
        sys.stderr.write(usage+"\n")
        exit(1)
    vcfilename=args[0]
    #maf=float(args[0])

    freqfh=open('freq.log', 'w')

    vcfh=open(vcfilename,'r')

    #instantiate a VcfFile object
    vcfobj=VcfFile(vcfilename)
    #parse its metainfo lines (ones that begin with ##)
    if options.noheader == False:
        vcfobj.parseMetaLines(vcfh)
    #vcfobj.printMetaLines()
    descriptors = vcfobj.getMetaInfoDescription()
    infoids=[]
    for (tag, description) in descriptors:
        infoids.append(tag)

    if options.maftag  not in infoids and options.maftag != 'QUAL' and options.noheader == False:
        sys.stderr.write(options.maftag + " tag not in ##INFO headers!\n")
        exit(1)

    if options.vtag  not in infoids and options.vtag != 'QUAL' and options.noheader==False:
        sys.stderr.write(options.vtag + " tag not in ##INFO headers!\n")
        exit(1)

   
    #vcfh.seek(0)
    if options.noheader == False:
        vcfobj.parseHeaderLine(vcfh)
  


    if options.variantype==None:
        variantpattern=options.vtag+'=(\w+);'
    else:
        variantpattern=options.vtag+'=('+options.variantype+');'
    mafpattern=options.maftag+'=(0.\d+)'

    #print mafpattern, variantpattern


    for dataline in vcfobj.yieldVcfDataLine(vcfh):
        #print dataline
        fields=dataline.strip().split('\t')

        (chrom,pos,id,ref,alt,qual,filtercode,info)=fields[0:8]
        #if filtercode != options.filter and options.filter != None : continue

        
        if re.search(variantpattern, info ) == None:
            #sys.stderr.write("no variant pattern\n")
            continue
        
        variant_type=re.search(variantpattern, info ).groups()[0]
        
        
        if re.search(mafpattern, info ) == None:
            #sys.stderr.write("No mafpattern!\n")
            #sys.stderr.write(dataline+"\n")
            continue
        
        maf_value=re.search(mafpattern, info ).groups()[0]
        
        if float(maf_value) <= options.leq and float(maf_value) >= options.geq:

            if options.quiet == False:
                print dataline
            logstring="\t".join([chrom,pos,id,ref,alt,variant_type, options.maftag, maf_value])
            freqfh.write(logstring+'\n')
        
        

if __name__ == "__main__":
    main()
