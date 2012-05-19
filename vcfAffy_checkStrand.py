#!/usr/bin/env python
import sys
import os
import string
import re
from optparse import OptionParser

import MySQLdb as mdb
from VcfFile import *
from common import *


def main():
    usage = "usage: %prog [options] file.vcf"
    parser = OptionParser(usage)
    parser.add_option("--table", type="string", dest="table", default="snpArrayAffy6", help="mysql table to select rows from default snpArrayAffy6 ")
    parser.add_option("--db", type="string", dest="db", default="hg19", help="db to use (default hg19)")
    parser.add_option("--server", type="string", dest="server", default="genome-mysql.cse.ucsc.edu")
    parser.add_option("--user", type="string", dest="user", default="genomep")
    parser.add_option("--passwd", type="string", dest="passwd", default="password")
    parser.add_option("--chr", type="string", dest="chrom", default=None)
    (options, args)=parser.parse_args()

    logfh=open('checkstrand.log', 'w')
    con = mdb.connect(options.server, options.user, options.passwd, options.db)

    vcfile=args[0]
    vcfh=open(vcfile, 'r')
    vcfobj=VcfFile(vcfh)
    vcfobj.parseMetaAndHeaderLines(vcfh)
    #print vcfobj.printMetaLines()

    for vrec in vcfobj.yieldVcfRecordwithGenotypes(vcfh):
        if options.chrom != None and vrec.getChrom() != options.chrom:
            pass
        rsid=vrec.getId()
    

        #rsid="rs3748597"
        with con:
            cur = con.cursor()
            cur.execute("select strand from " + options.table + "  where  rsId='"+rsid+"'")
            rows=cur.fetchall()
            strand=rows[0][0]
            #print strand, vrec.getId()

            if strand == '-':
                refComp=reversecomplement( vrec.getRef() )
                altComp=reversecomplement( vrec.getAlt() )
                logstring = "\t".join([rsid, vrec.getRef(), vrec.getAlt(), 'refComp',refComp, 'altComp',altComp])
                vrec.setRef(refComp)
                vrec.setAlt(altComp)
                logfh.write(logstring +"\n")
        print vrec.toStringwithGenotypes()

if __name__ == "__main__":
    main()
