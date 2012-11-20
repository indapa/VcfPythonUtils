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
    parser.add_option("--max", type="int", dest="max", help="skip records that are greater than or equal to max (default sys.maxint)", default=sys.maxint)
    #parser.add_option("--v", action="store_true", dest="snp",  help="restrict analysis to SNPs (must have INFO ID SNP in header")

    (options, args)=parser.parse_args()

    vcfilename=args[0]
    fileName, fileExtension = os.path.splitext(vcfilename)
    #nuller.12:80717441..80717681.vcf
    regionpattern='nuller.(\d+):(\d+)..(\d+)'
    results=re.search(regionpattern,fileName ).groups()
    regionstr="\t".join(list(results))
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
            if int(dp) >= options.max: continue
            depth_list.append(int(dp))

    maxDP=max( array (depth_list))
    minDP= min (array (depth_list))
    medianDP=median (array (depth_list))
    meanDP=mean( array(depth_list))
    length=len(depth_list)

    outstr="\t".join([regionstr, str(maxDP), str(minDP), str(medianDP), str(meanDP), str(length)])
    print outstr


if __name__ == "__main__":
    main()
