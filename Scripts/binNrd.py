#!/usr/bin/python
import sys
import os
import string
import re
import numpy as np
from optparse import OptionParser



def main():
    usage = "usage: %prog [options] arg"
    parser = OptionParser(usage)
    (options, args)=parser.parse_args()

    file=args[0]
    fh=open(file, 'r')

    headerline=fh.readline().strip()
    total_indiv=0

    nrd_omnimaf=[] # list in the form [ (nrd, omni_vcf) ,,, ] for each site

    bins = np.linspace(0, 1, 1000) #1000 bins between 0 and 1
    L =[ [ ] for i in range(0, len(bins)) ] # L is a list of lists [ [] ,,, [] ] 

    #print len(L), len(bins)

    #print bins
    
    
    


    #now iterate thru the bins, and place the site_nrd  in the appropriate lists in L
    # based on its omni_maf
    sys.stderr.write("binning site nrd based on omni maf ...\n")
    count=0
    for line in fh:
        #print line.strip()
        fields=line.strip().split('\t')

        omni_maf= float(fields[12])
        if fields[16] == 'NA': 
            continue
        site_nrd= float(fields[16])
        count+=1
        total_indiv = int(fields[5] )
        #print total_indiv, site_nrd, omni_maf, fields
        bin_number =-1
        for i in range(0, len(bins)-1 ):
            if omni_maf >= bins[i] and omni_maf < bins[i+1]:
                L[i].append(  site_nrd )
                #bin_number=i
                continue
        #if bin_number == -1: print bin_number, omni_maf, site_nrd
    sum=0 
    for i in range(0, len(L) ):
        if len(L) ==0: continue
        
        a= np.array( L[i])
        if len(a) == 0:
            continue
        print bins[i], np.median(a)

if __name__ == "__main__":
    main()
