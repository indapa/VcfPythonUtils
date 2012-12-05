#!/usr/bin/env python
import itertools
import gzip
import tabix
from intervalTree import *
import sys
from optparse import OptionParser

def returnIntervals(bedfh,size=10000000,overlap=0):
    regions=[]
    region_size=size
    for line in bedfh:
        if '_' in line: continue
        fields = line.strip().split("\t")

        chrom_name = fields[0]
        
        chrom_length = int(fields[2])
        region_start = 0
        
        while region_start < chrom_length-overlap:
            start = region_start
            end = region_start + region_size
            if end > chrom_length:
                end = chrom_length
            region_string = chrom_name + "\t" + str(region_start) + "\t" + str(end)
            regions.append(region_string)
            yield region_string
            region_start = end - overlap





def main():
    usage = "usage: %prog [options] query.vcf"
    parser = OptionParser(usage)
    #parser.add_option("--query",type="string", dest="query", help="Vcf file you want to query against a variant catalog", default=None)
    parser.add_option("--catalog", type="string", dest="catalog", help="Vcf file representing variant catalog (e.g. 1000G, HapMap, NHBLI-ESP)", default=None)
    parser.add_option("--chrombed", type="string", dest="chrombed", help="bed3 file with chromosome sizes", default="chromInfo.bed")
    (options, args)=parser.parse_args()
    print args
    if args[0] == None:
        sys.stderr.write("provide a query vcf file!\n")
        sys.exit(1)
    queryvcf=args[0]
    if options.catalog == None:
        sys.stderr.write("please provide a Vcf file for --catalog  option!\n")
        sys.exit(1)
    
    #tabix objects for query and catalog vcf
    tb_cat = tabix.Tabix(options.catalog)
    tb_query=tabix.Tabix(queryvcf)
    chroms=[]
    chrom_starts=[]
    chrom_ends=[]
    bedfh=open( options.chrombed, 'r')

    for region_string  in returnIntervals (bedfh):
        (chr, start, end)=region_string.strip().split('\t')
        print chr, str(start), str(end)
        sys.stderr.write("tabix processing for chromosome " + chr + "\n")

        catalog_features=[]

        query_results = [ fields for fields in tb_query.query(chr, 1, int(end)) ]
        catalog_results = [ fields for fields in tb_cat.query(chr, 1, int(end)) ]
    
        #make an interval tree for the catalog results?
        sys.stderr.write("collecting variant catalog features ...\n")
        for q in catalog_results:

            (chr, start, rsid, ref, alt)=q[0:5]
            catalog_features.append ( [chr, int(start)-1, int(start), rsid] )
        
        """ make an intervalTree object for the variants in the catalog """
        catalogTree= intervalTree( catalog_features, 1,2,0,int(end) )
        sys.stderr.write("querying intervalTree for overlaps with query file Vcf records ...\n")
        for q in query_results:
            (chr, start, rsid, ref, alt)=q[0:5]
            start=int(start)-1
            end=int(start)
            results = catalogTree.findRange([start,end])
        
            if len(results)  == 0: 
                print "\t".join(q)
        


if __name__ == "__main__":
    main()

