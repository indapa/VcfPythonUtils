#!/usr/bin/env python
import sys
import os
import string
import re
from optparse import OptionParser
from VcfFile import *
from numpy import *


""" given a VCF file determine the breadth vs depth. For each sample print a file with two columns of min depth & #bases with that min depth
    the number bins is set by option --maxdp (default is 40x)"""


def main():
    usage = "usage: %prog [options] file.vcf"
    parser = OptionParser(usage)
    parser.add_option("--info", type="string", dest="infotag", help="INFO tag site/sample depth", default="DP")
    parser.add_option("--filter", type="string", dest="filter", help="only analyze records with matching filter (default is PASS)", default="PASS")
    parser.add_option("--maxdp", type="int", dest="maxdp", default=40, help="max value of depth bin (default is 40x")
    (options, args)=parser.parse_args()


    breadth_array = zeros( (1,options.maxdp+1) )[0]

    vcfilename=args[0]
    vcfh=open(vcfilename,'r')

    #instantiate a VcfFile object
    vcfobj=VcfFile(vcfilename)
    #parse its metainfo lines (ones that begin with ##)
    vcfobj.parseMetaLines(vcfh)
    descriptors = vcfobj.getMetaFormatDescription()
    infoids=[]
    for (tag, description) in descriptors:
        infoids.append(tag)
    
    if options.infotag  not in infoids and options.infotag != 'QUAL':
        sys.stderr.write(options.infotag + " tag not in ##FORMAT headers!\n")
        exit(1)


    vcfh.seek(0)
    vcfobj.parseHeaderLine(vcfh)
    samplelist=vcfobj.getSampleList()

    breadth_array = zeros( (len(samplelist),options.maxdp+1) )
    print breadth_array

    print samplelist
    pattern=options.infotag+'=(\w+)'
    
    for vrec in vcfobj.yieldVcfRecordwithGenotypes(vcfh):
        if vrec.getFilter() != options.filter: continue
        print vrec.toString()
        zipgenolist=vrec.zipGenotypes(samplelist)
        for sample_index in range(0, len(zipgenolist)):
            (sample, geno) = zipgenolist[sample_index]
            sample_depth=0
            
            if geno.getFormatVal(options.infotag) == None:
                continue
            elif geno.getFormatVal(options.infotag) == '.':
                sample_depth=0
            else:
                sample_depth = int( geno.getFormatVal(options.infotag) )

            for depth_index in range(0,options.maxdp +1):
                if int(sample_depth) > depth_index:
                    breadth_array[sample_index][depth_index]+=1

    for i in range(0, len(samplelist)):
        print samplelist[i], len(breadth_array[i])
        dpfilename=samplelist[i]+".breadth.depth.dat"
        fh=open(dpfilename, 'w')
        for dp in range(0,options.maxdp+1):
           
            outstring= "\t".join([ str(dp),str(breadth_array[i][dp]) ])
            fh.write( outstring + "\n")
        
        

if __name__ == "__main__":
    main()
