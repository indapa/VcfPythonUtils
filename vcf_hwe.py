#!/usr/bin/python
import sys
import os
import string
import re
from optparse import OptionParser
from vcfIO import *

def main():
    """ given a VCF file compute the likelihood ratio test statistic for Hardy-Weinberg equilribum
        2 * sum obs_i * log(obs_i/exp_i)                                                        """

    usage = "usage: %prog [options] file.vcf"
    parser = OptionParser(usage)
    #parser.add_option("--name", type="string|int|boolean", dest="name", help="help string")

    (options, args)=parser.parse_args()

    vcf1_fname=args[0]
    vcf1_fh= open(vcf1_fname, 'r')
    vcf1_samples=get_vcfsamples(vcf1_fh)


    #reset the filehandle positions
    vcf1_fh.seek(0)

    header="\t".join(['chrom', 'position', 'lrts', 'obsAA', 'obsAB', 'obsBB', 'expAA', 'expAB', 'expBB'])
    print header

    for vcf1_line in get_vcfdataline_passfilter(vcf1_fh):
        vcf1_data=split_vcfdataline(vcf1_line)
        vcf1_formatstr= vcf1_data[8]
        vcf1_infostr = vcf1_data[7]

        vcf1_ziptuple=zip(vcf1_samples, vcf1_data[9::] )
        obs = return_observed_genotype_counts (vcf1_ziptuple)
        exp = return_expected_genotype_counts(vcf1_ziptuple)
        lrts=  hweLRT ( vcf1_ziptuple )

        print vcf1_data[0], vcf1_data[1], lrts, obs[0], obs[1], obs[2], exp[0], exp[1], exp[2]

if __name__ == "__main__":
    main()
