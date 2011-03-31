#!/usr/bin/python
import sys
import os
import string
import re
import numpy as np
from optparse import OptionParser



def main():
    usage = "usage: %prog [options] site.nrd.txt"
    parser = OptionParser(usage)
    (options, args)=parser.parse_args()

    file=args[0]
    fh=open(file, 'r')

    headerline=fh.readline().strip()
    total_indiv=0

    het_nrd=[]
    refhomoz_nrd=[]
    nonrefhomoz_nrd=[]

    for line in fh:
        fields=line.strip().split('\t')
        omni_maf= float(fields[12])

        if fields[16] !='NA': 
            het_nrd.append( float( fields[16] ) )

        if fields[15] !='NA': 
            refhomoz_nrd.append( float( fields[15] ) )

        if fields[17] !='NA': 
            nonrefhomoz_nrd.append( float( fields[17] ) )

        
    
    mean_het_nrd= np.mean ( np.array ( het_nrd ) )
    median_het_nrd= np.median ( np.array ( het_nrd ) )

    mean_homozref_nrd= np.mean ( np.array ( refhomoz_nrd ) )
    median_homozref_nrd= np.median ( np.array ( refhomoz_nrd ) )

    mean_homoznonref_nrd=  np.mean ( np.array ( nonrefhomoz_nrd ) )
    median_homoznonref_nrd= np.median ( np.array ( nonrefhomoz_nrd ) )


    print "\nmean homozygous reference discrepancy: ", str(mean_homozref_nrd)
    print "median homozygous reference discrepancy: ", str(median_homozref_nrd),"\n"

    print "mean homozygous non-reference discrepancy: ", str(mean_homoznonref_nrd)
    print "median homozygous non-reference discrepancy: ", str(median_homoznonref_nrd), "\n"

    print "mean het reference discrepancy: ", str(mean_het_nrd)
    print "median het  reference discrepancy: ", str(median_het_nrd),"\n"

if __name__ == "__main__":
    main()
