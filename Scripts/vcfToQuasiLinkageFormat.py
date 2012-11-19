#!/usr/bin/python
import sys
import os
import string
import re
from optparse import OptionParser
from VcfFile import *
from PedFile import *
from common import *

""" given a vcf file and a ped file print the GL to *.ped file with GLs for the genotypes """

def main():
    usage = "usage: %prog [options] file.vcf"
    parser = OptionParser(usage)
    parser.add_option("--ped", type="string", dest="pedfile", help="help string")
    
    (options, args)=parser.parse_args()
    mapfile=open('vcf.map', 'w')
    pedobj=Pedigree(options.pedfile)
    mapfh=open('vcf.map', 'w')
    vcfilename=args[0]
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
        chrom=int(vrec.getChrom())
        pos=int(vrec.getPos())
       
         
        id=vrec.getId()
        mapstr="\t".join([str(chrom), id, '0', str(pos)])
        mapfh.write(mapstr+"\n")
        genotypes=vrec.getGenotypes()
        g=[(0,0), (0,1), (0,2), (0,3), (1,1), (1,2), (1,3), (2,2), (2,3), (3,3)]
        gl=sorted ( [(order(x[0], x[1]), x) for x in g] )
        mapstr= "\t".join([vrec.getChrom(), vrec.getPos(), vrec.getRef(), vrec.getAlt()])
        mapfh.write(mapstr+"\n")
        (likelihoods)=map(lambda x: x.getFormatVal('GL'), genotypes)
        #print likelihoods, len(likelihoods)
        gls=zip(samplelist, likelihoods)
        
        
        #print gls
        for (sample, gl) in gls:
            genotypeSampleDict[sample].append(gl)
    
    for s in genotypeSampleDict.keys():
        #print s, len(genotypeSampleDict[s])
        pedobj.addGenotypes(genotypeSampleDict[s], s)
    pedobj.dumpToFile()
if __name__ == "__main__":
    main()
