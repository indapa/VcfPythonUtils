#!/usr/bin/env python
import sys
import os
import string
import re
from optparse import OptionParser

import math
""" just a simple script to help convert phred scaled probablity and error probablilty """
def main():
    usage = "usage: %prog [options] number"
    parser = OptionParser(usage)
    parser.add_option("--Q",  action="store_true",  dest="qscore", help="report the phred scaled quality")
    parser.add_option("--P", action="store_true", dest="pscore", help="report the error probability ")

    (options, args)=parser.parse_args()

    number=args[0]
    number=float(number)

    if options.pscore==True:
         print "error probability: ",   pow(10,(-number/10))

    elif options.qscore==True:
        print "phred-scaled probablity: ",   -10 * math.log10(number)

    else:
        sys.stderr.write("did you forget to specity --Q or --P option?\n")
        exit(1)
    
    


if __name__ == "__main__":
    main()
