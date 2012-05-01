#!/usr/bin/env python
from optparse import OptionParser
def returnIntervals(bedfh):
    regions=[]
    region_size=10000000
    for line in bedfh:
        if '_' in line: continue
        fields = line.strip().split("\t")
        chrom_name = fields[0]
        
        chrom_length = int(fields[2])
        region_start = 1
        
        while region_start < chrom_length:
            start = region_start
            end = region_start + region_size
            if end > chrom_length:
                end = chrom_length
            region_string = chrom_name + ":" + str(region_start) + ".." + str(end)
            regions.append(region_string)

            region_start = end

    return regions

def main():
    #freebayes.20120429
    usage = "usage: %prog [options] intervals.bed"
    parser = OptionParser(usage)
    parser.add_option("--prefix", type="string", dest="prefix", help="prefix of vcf files", default="freebayes.20120429")
    parser.add_option("--suffix", type="string", dest="suffix", help="suffix of vcf files", default=".vcf")

    (options, args)=parser.parse_args()
    outfile=options.prefix+"."+options.suffix


    bedfile=args[0]
    bedfh=open(bedfile, 'r')

    regions=returnIntervals(bedfh)
    
    files=[]

    for region in regions:
        region='chr'+region
        #print region
        vcf=".".join([options.prefix, region, options.suffix])
        files.append(vcf)

    listring=" --in ".join(files)
    print "vcf_concat.py  --in " + listring + " > " + outfile
if __name__ == "__main__":
    main()
