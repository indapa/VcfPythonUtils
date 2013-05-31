#!/usr/bin/env python
from itertools import *
from VcfFile import *
from VcfSampleEval import *
import numpy as np
import re
from optparse import OptionParser
import os
import pysam

def main():
    
    """ given a VCF file and bam file containing the sample(s) in the VCF this willl print out 
    a pileup count of the ref and alt allele that is in the VCF file """


    usage = "usage: %prog [option] file.vcf"
    parser =OptionParser(usage)
    parser.add_option("--bam", type="string", dest="bam", default=None, help="bam file to perform pileup on")
    parser.add_option("--mapq", type="float", dest="mapq", default=0., help="Exclude alignments from analysis if they have a mapping less than mapq (default is 0)")
    parser.add_option("--bq", type ="float", dest="bq", default =0. , help="Exclude bases from analysis if their supporting base quality is less that --bq (default is 0)")
    parser.add_option("--includeDuplicates", action="store_false", dest="duplicate", help="include duplicate marked reads in analysis (turned off by default) ")
    (options, args)=parser.parse_args()
    if options.bam == None:
        sys.stderr.write("please provide a value to --bam option\n")
        sys.exit(1)
    
    vcfilename=args[0]
    basename=os.path.splitext(vcfilename)[0]
    bamfilename=options.bam
    
    if os.path.exists(bamfilename+".bai") == False:
        sys.stderr.write("please check for existence of bam index file (*.bai)\n")
        exit(1)
        
    vcfobj=VcfFile(vcfilename)
    vcfh=open(vcfilename,'r')

    vcfobj.parseMetaAndHeaderLines(vcfh)
    header=vcfobj.returnHeader() +"\n"
        
    pybamfile = pysam.Samfile(bamfilename, "rb" )
    
    samples=vcfobj.getSampleList()
    
    print samples
    
    for vrec in vcfobj.yieldVcfRecord(vcfh):
        (chrom, start, end)=vrec.getChrom(), int( vrec.getPos() )-1, int(vrec.getPos() )
        print chrom, str(start), str(end)
        print vrec.getRef()
        
        for pileupcolumn in pybamfile.pileup( chrom, start, end):
            if pileupcolumn.pos != end:
                continue
            sys.stdout.write('chr'+chrom+ " " + str(start) +  " " + str(end) + " " + str(pileupcolumn.pos) + " ")
            print 'coverage at base %s = %s' % (pileupcolumn.pos , pileupcolumn.n)
            
            seqdict={}
            for (base,count) in ( ('A',0), ('C',0), ('G',0), ('T',0), ('N',0) ):
                seqdict[base]=count
            
            for pileupread in pileupcolumn.pileups:
                
                if pileupread.alignment.is_duplicate == True and options.duplicate == False: continue
                if pileupread.alignment.mapq < options.mapq: continue
                if  ( ord ( pileupread.alignment.qual[ pileupread.qpos -1 ] )  - 33 ) < options.bq: continue
                seqdict[ pileupread.alignment.seq[pileupread.qpos-1] ] +=1
                #print pileupread.alignment.seq, len(pileupread.alignment.seq), pileupread.qpos
            print vrec.getRef(), seqdict[vrec.getRef()]
            print vrec.getAlt(),seqdict[vrec.getAlt()]
            for nt in ('A', 'C', 'G', 'T', 'N'):
                sys.stdout.write( str(seqdict[nt]) + " ")
            sys.stdout.write("\n")
            
    pybamfile.close()
        
                
            
    
    
    
    
    
    
    


if __name__ == "__main__":
    main()

