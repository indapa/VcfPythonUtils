#!/usr/bin/python2.6
import sys
import os
import string
import re
from optparse import OptionParser
from VcfFile import *
from PedFile import *
from common import *

""" given a vcf file and a ped file print the genotypes to *.ped file and positions to a  corresponding map file """

def main():
    usage = "usage: %prog [options] file.vcf \n \ngenerate ped and map file from vcf and ped file with genotypes"
    parser = OptionParser(usage)
    parser.add_option("--ped", type="string", dest="pedfile", help="help string")
    parser.add_option("--filter", type="string", dest="filter", help="only analyze records with matching filter (default is None)", default=None)
    (options, args)=parser.parse_args()
   
    pedobj=Pedigree(options.pedfile)
   
    vcfilename=args[0]
    (name,extension)=os.path.splitext(vcfilename)
    mapfilename=name+".map"
    pedfilename=name+".ped"
    mapfh=open(mapfilename, 'w')
    pedfh=open(pedfilename, 'w')
    vcfh=open(vcfilename,'r')
    vcfobj=VcfFile(vcfilename)
    vcfobj.parseMetaLines(vcfh)
    vcfh.seek(0)
    
    vcfobj.parseHeaderLine(vcfh)
    genotypeSampleDict={}
    samplelist=vcfobj.getSampleList()
    for s in samplelist:
        genotypeSampleDict[s]=[]
    for vrec in  vcfobj.yieldVcfRecordwithGenotypes(vcfh):
        #print vrec.toString()
        filtercode = vrec.getFilter()
        if filtercode != options.filter and options.filter != None : continue
        chrom=int(vrec.getChrom())
        pos=int(vrec.getPos())
       
         
        id=vrec.getId()
        mapstr="\t".join([str(chrom), id, '0', str(pos)])
        mapfh.write(mapstr+"\n")
        genotypes=vrec.getGenotypes()
        #print genotypes
        genotype_strings=map(lambda x: x.getFormatVal('GT'), genotypes)
        
     
        gls=zip(samplelist, genotype_strings)
  
        for (sample, genos) in gls:
            genotypeSampleDict[sample].append(genos)
    
    for s in genotypeSampleDict.keys():
        pedobj.addGenotypes(genotypeSampleDict[s], s)

    pedobj.dumpToFile(pedfh)
if __name__ == "__main__":
    main()
