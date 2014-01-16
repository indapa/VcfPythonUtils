#!/usr/bin/env python
import sys
import os
import argparse
from PedPy import Ped, Pedfile
import gzip
from itertools import *
from VcfFile import *
from VcfSampleEval import *
from optparse import OptionParser
import argparse
import os



""" Given a gzipped vcf file and pedigree file, generate a new vcf with only those samples present in the pedigree (ped file)  """

def main():
    usage = "usage: %prog [options]  "
    parser = argparse.ArgumentParser(description='Given a gzipped vcf file and pedigree file, generate a new vcf with only those samples present in the pedigree (ped file) ')
    parser.add_argument('-ped', dest='pedfile', type=str, help="*.ped file")
    parser.add_argument('vcfile',  type=str,help='*.vcf.gz file')

    args=parser.parse_args()

    """ parse the pedfile and return the list of iids to keep from the VCF file """
    pedobj=Pedfile(args.pedfile)
    pedobj.parsePedfile()

    keeplist=  pedobj.returnIndivids()

    #open the VCFfile
    vcfh=gzip.open(args.vcfile,'r')
    vcfobj=VcfFile(args.vcfile)

    vcfobj.parseMetaAndHeaderLines(vcfh)
    samples=vcfobj.getSampleList()
    newsamples= [ s for s in samples if s in keeplist]

    print newsamples

    vcfobj.setSampleList(newsamples)
    header=vcfobj.returnHeader()
    print header

    for vrec in vcfobj.yieldVcfRecordwithGenotypes(vcfh):
        keepGenotypes=[]
        vrec_ziptuple=vrec.zipGenotypes(samples)
        for (s, genObj) in vrec_ziptuple:
            if s in keeplist:
                keepGenotypes.append( genObj )
    
        vrec.addGenotypeList(  keepGenotypes )
        print vrec.toStringwithGenotypes()



if __name__ == "__main__":
    main()
