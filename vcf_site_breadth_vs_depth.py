#!/usr/bin/env python
import sys
import os
import string
import re
from optparse import OptionParser
from VcfFile import *
from numpy import *

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
    descriptors = vcfobj.getMetaInfoDescription()
    infoids=[]
    for (tag, description) in descriptors:
        infoids.append(tag)

    if options.infotag  not in infoids and options.infotag != 'QUAL':
        sys.stderr.write(options.infotag + " tag not in ##INFO headers!\n")
        exit(1)


    vcfh.seek(0)
    vcfobj.parseHeaderLine(vcfh)

    pattern=options.infotag+'=(\w+)'

    for vrec in vcfobj.yieldVcfRecord(vcfh):
        if vrec.getFilter() != options.filter: continue
        if re.search(pattern, vrec.getInfo() ) == None: #regex search on the INFO column of the VCF datline record
            sys.stderr.write("regex to find site DP value failed!\n")
            print vrec.toString()
            exit(1)
        else:
            value=re.search(pattern, vrec.getInfo() ).groups()[0]
            #print value
            for i in range(0,options.maxdp +1):
                if int(value) > i:
                    breadth_array[i]+=1
        
    for i in range(0, options.maxdp+1):
        print i, breadth_array[i]

if __name__ == "__main__":
    main()
