#!/usr/bin/env python
import sys
import os
import string
import re
from optparse import OptionParser
from VcfFile import *
from Pedfile import *


def makeTupleStringRep(t):
    """http://stackoverflow.com/questions/1455602/printing-tuple-with-string-formatting-in-python"""
    s = "%s" % (t,)
    return s

def marginal_segregation( gtuple ):
    gamete_list=list(gtuple)
    if '.' in gamete_list:
        called_alleles =[ elem for elem in gamete_list if elem != '.' ]
        if len(called_alleles)==0:
            called_alleles=['0','1']
        for a in called_alleles:
            for i in ['0', '1']:
                yield ( a, i )


def doPunnett( maternal, paternal):
    m1,m2=maternal
    p1,p2=paternal
    genotype_space=[]
    for pallele in paternal:
        if pallele=='.': return None
        for mallele in maternal:
            if mallele == '.': return None
            genotype_space.append( pallele+mallele )
    return  list(set(genotype_space))

def doCrossAutosomal(maternal, paternal):
    if '.' in maternal and '.' not in paternal:
        for maternal_tuple in marginal_segregation(maternal):
            #print maternal_tuple, paternal
            yield doPunnett(maternal_tuple, paternal)
            #print "=="
    elif '.' in paternal and '.'  not in maternal:
        for paternal_tuple in marginal_segregation(paternal):
            yield doPunnett(paternal_tuple, maternal)
            #print paternal_tuple, maternal
            #print "=="
    elif '.' in paternal and '.' in maternal:
        for paternal_tuple in marginal_segregation(paternal):
            for maternal_tuple in marginal_segregation(maternal):
                yield doPunnett(paternal_tuple, maternal_tuple)
                #print paternal_tuple, maternal_tuple
                #print "=="
    else:
        yield doPunnett(maternal, paternal)
        #print paternal, maternal
        #print "=="



""" sets MENDEL for FILTER column  for Mendelian inconsistencies in genotypes of a VCF file for a nuclear family """
def main():
    usage = "usage: %prog [options] file.vcf\ncheck for Mendelian inconsistencies in genotypes of a VCF file\n"
    parser = OptionParser(usage)
    parser.add_option("--ped", type="string", dest="pedfile", help="ped format file of sample names that comprise the nucelar family", default="")
    parser.add_option("--info", type="string", dest="infotag", help="INFO tag id that annotates what type of variant the VCF record is", default="TYPE")
    parser.add_option("--type", type="string", dest="variantype", help="type of variant (SNP INS DEL) default is SNP", default="SNP")
    parser.add_option("--filter", type="string", dest="filter", help="only analyze VCF records matching FILTER (default is PASS)", default="PASS")
    (options, args)=parser.parse_args()

    if options.pedfile=="":
        sys.stderr.write("please provide a value for the --ped option!\n")
        exit(1)
    if args[0]==None:
        sys.stderr.write("please provide a vcf file!\n")
        exit(1)

    pattern=options.infotag+'=('+options.variantype+')'

    vcfobj=VcfFile(args[0])
    pedfileobj=Pedfile(options.pedfile)

    pedfh=open(options.pedfile, 'r')
    vcfh=open(args[0], 'r')
    mendelfh=open('mendel.log', 'w')

    pedfileobj.parsePedfile(pedfh)
    

    vcfobj.parseMetaLines(vcfh)
    vcfobj.addMetaFilterHeader("MENDEL", "mendelian inconsistency")
    vcfobj.printMetaLines()
    vcfh.seek(0)
    vcfobj.parseHeaderLine(vcfh)
    vcfobj.printHeaderLine()

    samplelist=vcfobj.getSampleList()
    founderlist=pedfileobj.returnFounderIds()
    nonfounderlist=pedfileobj.returnNonFounderIds()
    
    mendelfh.write("#nonfounder_sample chrom position nonfounder_alleles founder_alleles founder_alleles" +"\n")
    for vrec in vcfobj.yieldVcfRecordwithGenotypes(vcfh):

        searchresult=re.search(pattern, vrec.getInfo() )
        if re.search(pattern, vrec.getInfo() ) == None:
            print vrec.toStringwithGenotypes()
            continue

        if vrec.getFilter() != 'PASS':
            print vrec.toStringwithGenotypes()
            continue
        zipgenolist=vrec.zipGenotypes(samplelist)
             
        founderzipgenolist=[ (samp, gobj) for (samp, gobj) in zipgenolist if samp in founderlist ]
        nonfounderzipgenolist=[ (samp, gobj) for (samp, gobj) in zipgenolist if samp in nonfounderlist ]

        #print founderzipgenolist[0][0], founderzipgenolist[0][1].getAlleles()
        #print founderzipgenolist[1][0], founderzipgenolist[1][1].getAlleles()
        genotype_space=[]
        for possible_genotypes in doCrossAutosomal( founderzipgenolist[0][1].getAlleles(), founderzipgenolist[1][1].getAlleles()) :
            genotype_space+=possible_genotypes
        #print genotype_space
        for nonfounder,vcfgobj in nonfounderzipgenolist:
            nonfounder_gstring=vcfgobj.getAlleles()[0]+vcfgobj.getAlleles()[1]
             
            if nonfounder_gstring not in genotype_space and nonfounder_gstring[::-1] not in genotype_space and '.' not in nonfounder_gstring:
                outstring = "\t".join( [ nonfounder, vrec.getChrom(), vrec.getPos() , makeTupleStringRep ( vcfgobj.getAlleles() ), makeTupleStringRep( founderzipgenolist[0][1].getAlleles() ) , makeTupleStringRep (founderzipgenolist[1][1].getAlleles() ) ] )
                mendelfh.write( outstring +"\n")
                
                vrec.setFilter("MENDEL")
        print vrec.toStringwithGenotypes()
if __name__ == "__main__":
    main()
