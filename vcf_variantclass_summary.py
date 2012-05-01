#!/usr/bin/env python
import sys
import os
import string
import re
from optparse import OptionParser
from VcfFile import *

""" print the nuumber of each type of variant class ( e.g. snp insertion, deletion, mnp, complex in a VCF file """

def TsTvRatio ( vreclist ):
    ts=0
    tv=0
    
    for vrec in vreclist:
        if vrec.isTransition() != None:
            if vrec.isTransition() == True:
                ts+=1   
            else:
                tv +=1
    ratio = float(ts)/float(tv)
    return ratio


def main():
    usage = "usage: %prog [options] arg"
    parser = OptionParser(usage)
    parser.add_option("--info", type="string", dest="infotag", help="INFO tag id that annotates what type of variant the VCF record is", default="TYPE")
    parser.add_option("--filter", type="string", dest="filter", help="only analyze records with matching filter (default is None)", default=None)

    (options, args)=parser.parse_args()
    if options.infotag == "":
        sys.stderr.write("provide a value for --info parameter!\n")
        exit(1)


    variant_dict={} #key variant type value VcfRecord object

    vcfilename=args[0]
    vcfh=open(vcfilename,'r')

    #instantiate a VcfFile object
    vcfobj=VcfFile(vcfilename)
    #parse its metainfo lines (ones that begin with ##)
    vcfobj.parseMetaAndHeaderLines(vcfh)
    
    descriptors = vcfobj.getMetaInfoDescription()
    infoids=[]
    for (tag, description) in descriptors:
        tag
        infoids.append(tag)

    if options.infotag  not in infoids and options.infotag != 'QUAL':
        sys.stderr.write(options.infotag + " tag not in ##INFO headers!\n")
        exit(1)

    

    pattern=options.infotag+'=(\w+)'
    
    for vrec in vcfobj.yieldVcfRecord(vcfh):
        if vrec.getFilter() != options.filter and options.filter != None: continue
        
        searchresult=re.search(pattern, vrec.getInfo() )
        if re.search(pattern, vrec.getInfo() ) == None:
            continue
        else:
            value=re.search(pattern, vrec.getInfo() ).groups()[0]
            #rint value
            if value not in variant_dict.keys():
                variant_dict[value]=[]
                variant_dict[value].append( vrec )
            else:
                variant_dict[value].append( vrec )


    
    sum=0
    sys.stderr.write("types and count of different variant classes in " + vcfilename + "\n")
    for k in variant_dict.keys():
        print k, len( variant_dict[k] )
        sum+=len( variant_dict[k] )
    print "TOTAL:", sum
    

if __name__ == "__main__":
    main()
