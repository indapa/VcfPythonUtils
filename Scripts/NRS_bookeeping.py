#!/usr/bin/env python
import sys
import os
import string
import re
from optparse import OptionParser
from VcfFile import *
from collections import *

def main():

    usage = "usage: %prog [options] file.vcf\n"
    parser = OptionParser(usage)

    (options, args)=parser.parse_args()

    maxdel=0.05 #Maximum fraction of reads with deletions spanning this locus for it to be callable in GATK UF
    nrs_counter=Counter()

    vcfilename=args[0]
    vcfh=open(vcfilename, 'r')

    vcfobj=VcfFile(vcfilename)
    vcfobj.parseMetaAndHeaderLines(vcfh)

    dp_pattern='DP=(\d+)'
    mq_pattern='MQ0=(\d+)'
    del_pattern='Dels=(0.\d+)'


    for vrec in vcfobj.yieldVcfRecordwithGenotypes(vcfh):
        if len( vrec.getRef() ) != 1 or len( vrec.getAlt() ) !=1:
            nrs_counter[ 'indel' ] += 1
            continue
    
        if vrec.getFilter() != 'PASS':
            nrs_counter[ vrec.getFilter() ] += 1
            continue
    
        
    
        if re.search(dp_pattern, vrec.getInfo() ) == None:
            nrs_counter[ 'no read depth' ] +=1
            continue
        if re.search( del_pattern, vrec.getInfo() ) != None:
            del_frac=re.search(del_pattern, vrec.getInfo() ).groups()[0]
            if float(del_frac) >=maxdel:
                nrs_counter[ 'maxdel fraction']+= 1
                continue
            else:
            #print vrec.toStringwithGenotypes()     
                dp=re.search(dp_pattern, vrec.getInfo() ).groups()[0]
                mqzero=re.search(mq_pattern, vrec.getInfo() ).groups()[0]
            #print str(dp), str(mqzero)
            
                if float(mqzero)/float(dp) > .10:
                    nrs_counter[ 'mqzero' ] +=1
                    continue
                else:
                    nrs_counter['no read depth']+=1
                    print vrec.toStringwithGenotypes()
                    
            
                
        
    for (type,count) in nrs_counter.items():
        print type, count
    print sum(nrs_counter.values()) 


if __name__ == "__main__":
    main()


