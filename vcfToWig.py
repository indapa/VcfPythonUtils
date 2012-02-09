#!/usr/bin/python
import sys
import os
import string
import re
from optparse import OptionParser
from VcfFile import *
from numpy import *
import glob
""" convert vcf to wig format to plot a numerical attribute in UCSC browser using wig format"""
""" http://genome.ucsc.edu/goldenPath/help/wiggle.html  """
""" assumes vcf has data from a single chromosome on """

def main():
    usage = "usage: %prog [options] file.vcf"
    parser = OptionParser(usage)
    parser.add_option("--info", type="string", dest="infotag", help="INFO tag site/sample depth", default="DP")
    parser.add_option("--filter", type="string", dest="filter", help="only analyze records with matching filter (default is .)", default=".")
    
    parser.add_option("--chrom", type="string", dest="chr", default="chrN")
    parser.add_option("--override", action="store_true", dest="override", default=False,  help="delete all previously existing wig files in directory before writing new ones")
    (options, args)=parser.parse_args()
    print options
    vcfilename=args[0]
    vcfh=open(vcfilename,'r')
    fileName, fileExtension = os.path.splitext(vcfilename)
    sys.stderr.write("processing " + fileName + "\n")
   
    #instantiate a VcfFile object
    vcfobj=VcfFile(vcfilename)
    #parse its metainfo lines (ones that begin with ##)
    vcfobj.parseMetaLines(vcfh)
    descriptors = vcfobj.getMetaFormatDescription()
    infoids=[]
    for (tag, description) in descriptors:
        infoids.append(tag)

    if options.infotag not in infoids and options.infotag != 'QUAL':
        sys.stderr.write(options.infotag + " tag not in ##FORMAT headers!\n")
        exit(1)


    vcfh.seek(0)
    vcfobj.parseHeaderLine(vcfh)
    samplelist=vcfobj.getSampleList()
    wigs={}
    for s in samplelist:
        wigfile=s+".wig"
        if os.path.exists(wigfile):
            if not options.override:
                wigs[s]=open(wigfile, 'a')
            else:
                wigs[s]=open(wigfile, 'w')
                wigs[s].write("track type=wiggle_0 name=\"" + s + "\" visibility=full\n")
                wigs[s].write("variableStep  chrom="+options.chr+"\n")
        else:
            wigs[s]=open(wigfile, 'w')
            wigs[s].write("track type=wiggle_0 name=\"" + s + "\" visibility=full\n")
            wigs[s].write("variableStep  chrom="+options.chr+"\n")
    
    
    pattern=options.infotag+'=(\w+)'

    for vrec in vcfobj.yieldVcfRecordwithGenotypes(vcfh):
        if vrec.getFilter() != options.filter: continue
        zipgenolist=vrec.zipGenotypes(samplelist)
        depths=[]
        for (sample, geno) in zipgenolist:
            
            sample_depth =  geno.getFormatVal(options.infotag)
            if sample_depth == '.':
                sample_depth=0
            
            outstr=" ".join(  [ str( vrec.getPos() ), str(sample_depth)  ] )
            wigs[sample].write(outstr+"\n")
  

    for s in samplelist:
        wigs[s].close()
        

    
if __name__ == "__main__":
    main()
