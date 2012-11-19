#!/usr/bin/python

from optparse import OptionParser
from bx.intervals import *
import sys

def genotype_emmissions (ziptuple, intersected_ziptuple):
    a1, a2, a1_i, a2_i = (".", ".", ".", ".")
    emitted_gt=""
    emitted_tuple = []

    for i in range(0, len(ziptuple)):
        (sample, genotype) = ziptuple[i]
        (gt, dp) = genotype.split(':')

        (sample_i, genotype_i) = intersected_ziptuple[i]
        (gt_i, dp_i) = genotype_i.split(':')

        if gt == '.':
            (a1, a2) = ('.', '.')
        else:
            (a1, a2) = gt.split('/')

            if gt_i == '.':
                (a1_i, a2_i) = ('.', '.')
            else:
                (a1_i, a2_i) = gt_i.split('/')


                #print sample, gt, a1, a2, "\t",  gt_i, a1_i, a2_i

        if a1 == a2 and a1_i == a2_i: #both homoz
            #print "homoz"
            if a1 == a2  == a1_i == a2_i: #are the matching homoz?
                emitted_gt = ":".join([gt,dp])
            else:
                #print "homoz mismatch"
                emitted_gt = ":".join([".",dp])

            if a1 != a2 and a1_i != a2_i:
                #print "het"
                emitted_gt = ":".join([gt,dp])

            if a1 != a2 and a1_i == a2_i:
                #print "het/homz mismatch"
                emitted_gt = ":".join([".",dp])
            if a1 == a2 and a1_i != a2_i:
                #print "het/homz mismatch"
                emitted_gt = ":".join([".",dp])
        emitted_tuple.append( (sample, emitted_gt))

        emissions=[]
        for (sample, gt) in emitted_tuple:
            emissions.append(gt)
    return emissions

def makeBedInterval(bedfh):
    sys.stderr.write("making bed intervals ....\n")
    bed_ranges = {}
    
    for line in bedfh:
        if '#' in line:
            continue
        fields = line.strip().split('\t')
        chrom, start, end = fields[0], fields[1], fields[2]
       
        otherstring = "\t".join(fields[3:])
        if chrom not in bed_ranges: bed_ranges[chrom] = Intersecter()
        bed_ranges[chrom].add_interval(Interval(int(start), int(end), strand=1, info={'chr':chrom, 'other':otherstring}) )
    
    return bed_ranges

def makeVcfInterval(vcffh):

    sys.stderr.write("making vcf intervals ...\n")
    vcf_ranges = {}
    names=[]
    fh = open('vcf.bed', 'w+')
    total_vcfrecords =0
    for line in vcffh:
        if '#CHROM' in line:
            fields = line.strip().split('\t')
            names = fields[9::]

        if '#' in line:
            continue#print line.strip()

        fields = line.strip().split('\t')
        genotypes = fields[9::]
        chrom,start,id,ref,alt,qual,filter,info = fields[0], int(fields[1]), fields[2], fields[3], fields[4], fields[5], fields[6], fields[7]
        #print filter
        if filter != "0": continue
        #if filter != "PASS": continue
        
        chrom = "chr"+chrom
        otherstring = "\t".join(fields[9:])
        
        bedout = "\t".join([chrom, str(start-1), str(start), "\n"])
        fh.write(bedout)
        if chrom not in vcf_ranges: vcf_ranges[chrom] = Intersecter()
        vcf_ranges[chrom].add_interval( Interval( start-1, start,strand=1, name=id, info={'chr':chrom, 'ref':ref, 'alt':alt, 'qual':qual, 'filter':filter, 'info':info, 'samples':names, 'genotypes':genotypes} ) )
        total_vcfrecords = total_vcfrecords+1
    #print vcf_ranges
    return  vcf_ranges

def main():
    usage = "usage: %prog [options] arg"
    parser = OptionParser(usage)


    parser.add_option("--vcf", type = "string", dest= "vcfile", help="vcf file to analzye")
    parser.add_option("--bed", type = "string", dest = "bedfile", help ="bed file to analyze")
    parser.add_option("--chrom", type = "string", dest="chromosome", help="chrmosome to analyze")
    parser.add_option("--v", action="store_true", dest="diff", default = False, help="elements in the bed not in the vcf")

    (options, args) = parser.parse_args()
   
    
    if len(args) != 1:
        parser.error("incorrect number of arguments")
        parser.usage()

    infile = args[0]

    if options.bedfile != None:
        try:
            sys.stderr.write("opening bed file\n")
            bed_fh =open(options.bedfile)
        except:
            sys.stderr.write("cannot open bed file!\n")
            return
        interval_ranges = makeBedInterval(bed_fh)

    if options.vcfile != None:
        try:
            sys.stderr.write("opening vcf file \n")
            vcf_fh = open(options.vcfile)
        except:
            sys.stderr.write("cannot open vcf file\n")
            return
        interval_ranges = makeVcfInterval(vcf_fh)

    try:
        fh = open(infile)
    except:
        sys.stderr.write("cannot open infile!\n")
        return
    names = []
    sys.stderr.write("intersecting vcf with interval ranges ...\n")
    for line in fh:
        if '#CHROM' in line:
            fields = line.strip().split('\t')
            names = fields[9::]
        if '#' in line:
            print line.strip()
            continue
        fields = line.strip().split('\t')
        metadata = fields[0:9]
        chr, start = fields[0], int(fields[1])
        #chr = "chr"+chr
        if chr == "X": continue
        if chr == "23": continue
        genotypes = fields[9::]
        ziptuple = zip(names, genotypes)
        if chr not in interval_ranges.keys():
            if options.diff == True:
                print line.strip()
            else: continue

        if len( interval_ranges[chr].find( start-1, start ) ) > 0 and options.diff == False:
            print line.strip()
            continue
            #intersectobj = interval_ranges[chr].find( start-1, start )

            #intersected_genotypes = intersectobj[0].info['genotypes']
            #intersected_samples = intersectobj[0].info['samples']
            
            #intersected_ziptuple = zip(intersected_samples, intersected_genotypes)
            
            #metafields = "\t".join(fields[0:9])
            #outstr = "\t".join(emissions)
            #print metafields, "\t", outstr

        if len( interval_ranges[chr].find( start-1, start ) ) <= 0 and options.diff == True:
           print line.strip()
if __name__ == "__main__":
    main()
