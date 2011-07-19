#!/usr/bin/env python
import sys
import os
import string
import re
from optparse import OptionParser
from numpy import *
from VcfFile import *
import rpy2.robjects as robjects
from rpy2.robjects.packages import importr


""" print some summary statistics of depth at sites in a VCF file and generate an Rplot of the percentiles """

def main():
    usage = "usage: %prog [options] arg"
    parser = OptionParser(usage)
    usage = "usage: %prog [options] file.vcf \n print summary information about site depth in records of a VCF file\n"
    parser = OptionParser(usage)
    #parser.add_option("--v", action="store_true", dest="snp",  help="restrict analysis to SNPs (must have INFO ID SNP in header")

    (options, args)=parser.parse_args()

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

    if 'DP' not in infoids:
        sys.stderr.write("DP tag not in ##INFO headers!")
        exit(1)

    vcfh.seek(0)
    vcfobj.parseHeaderLine(vcfh)

    pattern='DP=(\d+)'
    depth_list=[]
    for vrec in vcfobj.yieldVcfRecord(vcfh):

        dp=re.search(pattern, vrec.getInfo() ).groups()[0]
        if dp == None:
            sys.stderr.write("unable to parse DP value from INFO field\n")
            continue
        else:
            depth_list.append(int(dp))
    print "max DP: ", max( array (depth_list))
    print "min DP: ", min (array (depth_list))
    print "median DP: ", median (array (depth_list))
    print "mean DP: ", mean( array(depth_list))
    print "total sites analyzed: " ,len(depth_list)


    grdevices = importr('grDevices')
    grdevices.png(file="depth.png",width=512, height=512)

    r = robjects.r
   

    x = robjects.IntVector(depth_list)
    ecdf = robjects.r['ecdf']
    Fn = ecdf(x)
    percentiles=Fn(x)
    
    r.plot(percentiles,x, xlab="percentile", ylab="depth", main="Depth Percentile")
    grdevices.dev_off()

    sys.stderr.write("see R plot in depth.png")

if __name__ == "__main__":
    main()
