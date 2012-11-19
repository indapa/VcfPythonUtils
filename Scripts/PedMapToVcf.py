#!/usr/bin/env python
import sys
import os
import string
import re
from optparse import OptionParser
from VcfFile import *
from PedFile import *
import numpy as np

""" convert to VCF file, given a ped/map file """
# bioinformatics is so much fun!
def main():
    usage = "usage: %prog [options] file\n\nconvert file.ped file.map to file.vcf"
    parser = OptionParser(usage)
    #parser.add_option("--name", type="string|int|boolean", dest="name", help="help string")

    (options, args)=parser.parse_args()

    filestub=args[0]
    pedfile=filestub+".ped"
    mapfile=filestub+".map"

    if os.path.exists(pedfile) == False:
        sys.stderr.write(pedfile + " doesn't exist!\n")
        exit(1)

    if os.path.exists(mapfile) == False:
        sys.stderr.write(mapfile + " doesn't exist!\n")
        exit(1)

    mapobjects=[]

    mapfh=open(mapfile, 'r')
    for line in mapfh:
        mapobj=Map(line.strip())
        mapobjects.append(mapobj)

    pedobj=Pedigree(pedfile)
    # initially the rows are the individuals, column are the markers
    # we transpose on the fly and make the rows the genotypes and the columns the individuals, since this is the structure of the VCF
    genotype_matrix=np.array( pedobj.getGenotypeMatrix() ).transpose()
    names=genotype_matrix[0,:]

   
    headerstring="\t".join(['#CHROM', 'POS', 'ID', 'REF', 'ID', 'QUAL', 'FILTER', 'INFO', 'FORMAT'])
    print headerstring+"\t" + "\t".join(names)
    (nrow, ncol) = np.shape(genotype_matrix)
    
 
    if len(mapobjects) != nrow-1:
        sys.stderr.write("unequal number of genotype/indiv and number of markers in map file!\n")
        exit(1)

    #now we print the #CHROM header
    # the trouble is information like the  in a ped/map file

    for i in range(nrow):
        if i==0: continue
        gls=zip(names, list(genotype_matrix[i,:]))
        gts=[]
        for (name, gt) in gls:
            gt=gt.replace('-1 -1', '.')
            gt=gt.replace(' ', '/')
            gts.append(gt)
        gstring="\t".join(gts)
        vcfcols="\t".join(['0', '1', '.', '.', '.', 'GT'])
        
        vcfstr=mapobjects[i-1].toStringVcf()+ "\t" + vcfcols + "\t" + gstring
        print vcfstr



if __name__ == "__main__":
    main()
