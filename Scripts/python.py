#!/usr/bin/python
import sys
import os
import string
import re
from optparse import OptionParser

def main():
    usage = "usage: %prog [options] arg"
    parser = OptionParser(usage)
    parser.add_option("--name", type="string|int|boolean", dest="name", help="help string")

    (options, args)=parser.parse_args()



if __name__ == "__main__":
    main()
