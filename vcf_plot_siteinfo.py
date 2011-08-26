#!/usr/bin/env python
import sys
import os
import string
import re
from optparse import OptionParser
from VcfFile import *
import rpy2.robjects as robjects
from rpy2.robjects.packages import importr


""" print 5 number summary (min, max, first & third quartile, mean, and median of   INFO (as well as QUAL)  tags  attribute at sites in a VCF file and generate an Rplot of the percentiles """

def main():
    usage = "usage: %prog [options] file.vcf"
    parser = OptionParser(usage)
    usage = "usage: %prog [options] file.vcf \n print 5 number summary  information about --info tag ids  in records of a VCF file\n"
    parser = OptionParser(usage)
    parser.add_option("--noplot", action="store_true", dest="noplot", help="only report 5 number summary; do not make R plot of percentiles")
    parser.add_option("--boxplot", action="store_true", dest="boxplot", help="make boxplot")
    parser.add_option("--dump", action="store_true", dest="dump", help="dump data to be plotted/summarized to file")
    parser.add_option("--info", type="string", dest="infotag", help="INFO tag attribute to get stats for", default="")
    parser.add_option("--max", type="int", dest="max", help="skip records that are greater than or equal to max (default sys.maxint)", default=sys.maxint)

    
    (options, args)=parser.parse_args()
    if options.infotag == "": 
        sys.stderr.write("provide a value for --info parameter!\n")
        exit(1)

    vcfilename=args[0]
    vcfh=open(vcfilename,'r')

    #instantiate a VcfFile object
    vcfobj=VcfFile(vcfilename)
    #parse its metainfo lines (ones that begin with ##)
    vcfobj.parseMetaLines(vcfh)
    descriptors = vcfobj.getMetaInfoDescription()
    infoids=[]
    for (tag, description) in descriptors:
        infoids.append(tag)

    if options.infotag  not in infoids and options.infotag != 'QUAL':
        sys.stderr.write(options.infotag + " tag not in ##INFO headers!\n")
        exit(1)

    vcfh.seek(0)
    vcfobj.parseHeaderLine(vcfh)

    pattern=options.infotag+'=(\d+?\.?\d+)'

    infovalues=[]
    for vrec in vcfobj.yieldVcfRecord(vcfh):
        

        if options.infotag == 'QUAL':
            if float(vrec.getQual() ) >= float(options.max):
                continue
            else:
                infovalues.append(float(vrec.getQual() ) )
        else:
            searchresult=re.search(pattern, vrec.getInfo() )
            if re.search(pattern, vrec.getInfo() ) == None:
                continue
            else:
                value=re.search(pattern, vrec.getInfo() ).groups()[0]
                if float(value) >= float(options.max): continue
                infovalues.append(float(value))

    
   


    grdevices = importr('grDevices')
    plotname=options.infotag+".png"
    grdevices.png(file=plotname,width=512, height=512)

    r = robjects.r
   

    x = robjects.FloatVector(infovalues)
    
    summary=robjects.r['summary']
    print summary(x)
    if options.noplot==True:
        return
    ecdf = robjects.r['ecdf']
    Fn = ecdf(x)
    percentiles=Fn(x)
    filename=vcfilename + "\n percentile plot " + options.infotag 
    sys.stderr.write("making R plot ...\n")
    r.plot(percentiles,x, xlab="percentile", ylab=options.infotag, main=filename)
    
    grdevices.dev_off()
    sys.stderr.write("see R plot in "+ plotname+"\n")
    if options.boxplot ==True:
        boxplotname=options.infotag+".boxplot.png"
        grdevices.png(file=boxplotname,width=512, height=512)
        r.boxplot(x, ylab=options.infotag, main=vcfilename + "\n boxplot " + options.infotag)

        grdevices.dev_off()
        sys.stderr.write("see R boxplot in "+ options.infotag+".boxplot.png"+"\n")


    if options.dump == True:

        dumpfile= vcfilename + "." + options.infotag + ".dat"
        dumpfh=open(dumpfile, 'w')
        sys.stderr.write("dumping data to " + dumpfile + "\n")
        for val in infovalues:
            dumpfh.write(str(val)+"\n")

if __name__ == "__main__":
    main()
