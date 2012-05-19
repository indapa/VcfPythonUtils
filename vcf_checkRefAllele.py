#!/usr/bin/env python
import sys
import os
import string
import re
from optparse import OptionParser
from VcfFile import *
import bx.seq.twobit

""" given a twobit file check to see if the ref allele matches the 2bit at the given position for snp variants
    if not then swap ref in vcf for the one from the 2bit"""
def main():
    usage = "usage: %prog [options] arg"
    parser = OptionParser(usage)
    parser.add_option("--tbf", type="string", dest="tbf", help="2bit file")
    parser.add_option("--info", type="string", dest="infotag", help="INFO tag id that annotates what type of variant the VCF record is", default="TYPE")
    parser.add_option("--type", type="string", dest="variantype", help="type of variant (SNP INS DEL)", default='snp')

    (options, args)=parser.parse_args()
    
    logfh=open('ref.log', 'w')
    
    try:
        sys.stderr.write("opening twobitfile...\n")
        twobit=bx.seq.twobit.TwoBitFile( open( options.tbf ) )
    except:
        sys.stderr.write("unable to open twobit file!\n")
        exit(1)

    #open vcf
    vcfile=args[0]
    vcfh=open(vcfile, 'r')
    vcfobj=VcfFile(vcfh)
    vcfobj.parseMetaAndHeaderLines(vcfh)
    pattern=options.infotag+'=('+options.variantype+')'
    
    for vrec in vcfobj.yieldVcfRecordwithGenotypes(vcfh):
        pos=vrec.getPos()
        start=int(pos)-1
        end=int(pos)
        
        info=vrec.getInfo()
        if re.search(pattern, info ) == None:
            continue
        else:
            value=re.search(pattern, info ).groups()[0]
            #print vrec.toString()
            assert( end > start ),"end greater than start!"
            try:
                sequence=twobit['chr'+vrec.getChrom()][start:end]
                sequence=sequence.upper()
            except:
                sys.stderr.write("unable to fetch sequence from 2bit file!\n")
        if sequence != vrec.getRef() and sequence == vrec.getAlt():
            vrec.setAlt( vrec.getRef() )
            vrec.setRef(sequence)
            logfh.write(vrec.getId() +"\n")
    print vrec.toStringwithGenotypes()
if __name__ == "__main__":
    main()
