#!/usr/bin/env python
import sys
import os
import string
import re
from optparse import OptionParser

from VcfFile import *
import bx.seq.twobit

""" given a 2bit file and a VCF extract out sequence from the 2bit file. Please note 2bit file extracton works on zero-based, half-open interval
    please note, as of now only SNPs are extracted. Complex variant types like indels, MNPs, SVs are ignored for now...."""

def main():
    usage = "usage: %prog [options] file.vcf"
    parser = OptionParser(usage)
    parser.add_option("--tbf", type="string", dest="tbf", help="2bit file")
    parser.add_option("--pad", type="int", dest="pad", default=0, help="extract sequence  upstream and downstream of position by pad value")

    parser.add_option("--info", type="string", dest="infotag", help="INFO tag id that annotates what type of variant the VCF record is", default="TYPE")
    parser.add_option("--type", type="string", dest="variantype", help="type of variant (SNP INS DEL)", default='snp')


    (options, args)=parser.parse_args()
    #open 2bitfile
    try:
        sys.stderr.write("opening twobitfile...\n")
        twobit=bx.seq.twobit.TwoBitFile( open( options.tbf ) )
    except:
        sys.stderr.write("unable to open twobit file!\n")
        exit(1)


    #open the vcf file
    vcfile=args[0]
    vcfh=open(vcfile, 'r')
    vcfobj=VcfFile(vcfh)
    vcfobj.parseMetaAndHeaderLines(vcfh)
    pattern=options.infotag+'=('+options.variantype+')'


    sequence=''
    downstream_seq=''
    upstream_seq=''

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
        
            if options.pad !=0:
                downstream_start=int(pos)
                upstream_end=int(pos)-1
                downstream_seq=twobit['chr'+vrec.getChrom()][downstream_start:downstream_start+options.pad]
                upstream_seq=twobit['chr'+vrec.getChrom()][upstream_end-options.pad:upstream_end]
                outstr="\t". join(['chr'+vrec.getChrom(), str(start), str(end), sequence, str(upstream_end-options.pad), str(upstream_end), upstream_seq, str(downstream_start), str(downstream_start+options.pad),downstream_seq] )
            else:
                outstr="\t". join(['chr'+vrec.getChrom(), str(start), str(end), sequence])
            print outstr

if __name__ == "__main__":
    main()
