#!/usr/bin/env python
import sys
import os
import string
import re
from optparse import OptionParser
from VcfFile import *
from PedFile import *
from common import *

""" given a vcf file generate a MendelSoft *.pre input file LINKAGE format described here http://www.inra.fr/mia/T/MendelSoft/ """

def main():
    usage = "usage: %prog [options] file.vcf \n \ngenerate Mendelsoft LINKAGE format input file"
    parser = OptionParser(usage)
    parser.add_option("--ped", type="string", dest="pedfile", help="help string")
    parser.add_option("--chrom", type="string", dest="chrom", default=None, help="restrict to chr")
    (options, args)=parser.parse_args()

    pedobj=Pedigree(options.pedfile)
    pediids=pedobj.getiIds()
    

    vcfilename=args[0]
    (name,extension)=os.path.splitext(vcfilename)
    mapfilename=name+".map"
    pedfilename=name+".pre"
    mapfh=open(mapfilename, 'w')
    pedfh=open(pedfilename, 'w')
    vcfh=open(vcfilename,'r')
    vcfobj=VcfFile(vcfilename)
    vcfobj.parseMetaLines(vcfh)
    vcfh.seek(0)

    vcfobj.parseHeaderLine(vcfh)
    genotypeSampleDict={}
    samplelist=vcfobj.getSampleList()
    samplesNotInVcf= [ s for s in pediids if s not in samplelist]
    
    
    for vrec in  vcfobj.yieldVcfRecordwithGenotypes(vcfh):
        #print vrec.toString()

        for s in pediids:
            genotypeSampleDict[s]=[]

        chrom=int(vrec.getChrom())
        if options.chrom != None and chrom != int(options.chrom):
            continue

        pos=int(vrec.getPos())


        id=vrec.getId()
        mapstr="\t".join([str(chrom), id, '0', str(pos)])
        prefilename=".".join(["chr",str(chrom), str(pos),"pre"])
        prefh=open(prefilename,'w')
        print prefilename
        mapfh.write(mapstr+"\n")
        genotypes=vrec.getGenotypes()
        #print genotypes
        #genotype_strings=map(lambda x: x.getFormatVal('GT'), genotypes)
        genotype_strings=map(lambda x: x.getAlleles(), genotypes)
        isCalled=map(lambda x: x.isCalled(), genotypes)
        #print isCalled
        if False in isCalled:
            sys.stderr.write("skipped site " + str(chrom) + " " + str(pos) + "\n")
            continue
            
        genotype_quals=map(lambda x: x.getFormatVal('GQ'), genotypes)
        #print genotype_quals
        #print genotype_strings
        genotype_ints= [map(int, x) for x in genotype_strings ]
        
        genotype_newstrings=[]
        for lst in genotype_ints:
            for i in range(len(lst)):
                lst[i]+=1
            
            genotype_newstrings.append( " ".join( [ str(x) for x in lst ] ) )
       
        

        gls=zip(samplelist, genotype_newstrings, genotype_quals)

        for (sample, genos, qual) in gls:
        #    genos=genos+' ' + qual
            genotypeSampleDict[sample].append(genos)
        for sample in samplesNotInVcf:
           genotypeSampleDict[sample].append('0 0')

        for s in genotypeSampleDict.keys():
            pedobj.addGenotypes(genotypeSampleDict[s], s)

        pedobj.dumpToFileNoPheno(prefh)
        print "=="
if __name__ == "__main__":
    main()

