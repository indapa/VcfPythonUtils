#!/usr/bin/env python
from itertools import *
from VcfFile import *
import numpy as np
import re
from optparse import OptionParser
from common import grouper
from common import typeofGenotype
import os

def main():
    usage = "usage: %prog [options] file.vcf \n calcuate NRS and NRD on a vcf generated from CombineVariants --genotypemergeoption UNIQUIFY\n"
    parser = OptionParser(usage)
    parser.add_option("--matrixonly", action="store_true", dest="matrixonly", help="only print concordance matrixe", default=False)
    parser.add_option("--includeRef", action="store_true", dest="includeRef", help="include sites in the set ReferenceInAll", default=False)

    (options, args)=parser.parse_args()


    vcfilename=args[0]
    basename=os.path.splitext(vcfilename)[0]
#row is veal column is comparison
    concordancetable= np.matrix( [ [ 0,0,0,0 ], [ 0,0,0,0 ], [ 0,0,0,0 ], [ 0,0,0,0 ] ] )
    calledtable = np.matrix ( [ [0 ,0] , [0,0] ] )
    #log file of sites that contribute to NRS penalty; hom-ref and no-calls at variant sites in comparison set
    nrslog=".".join([basename, 'nrs','log'])
    nrdlog=".".join([basename, 'nrd','log'])
    filterlog=".".join([basename, 'filtered','log'])
    multialleliclog=".".join([basename, 'multiallelic','log'])
    concordancelog=".".join([basename, 'concordance','log'])
    fieldslog=".".join([basename, 'fields', 'log'])
    nrsfh=open(nrslog, 'w')
    nrdfh=open(nrdlog, 'w')
    filteredfh=open(filterlog, 'w')
    multifh=open(multialleliclog, 'w')
    concordancefh=open(concordancelog, 'w')
    fieldsfh=open(fieldslog, 'w')
    fieldsfh.write('set'+"\n")
    vcfobj=VcfFile(vcfilename)
    vcfh=open(vcfilename,'r')

    vcfobj.parseMetaAndHeaderLines(vcfh)
    header=vcfobj.returnHeader() +"\n"
    nrsfh.write(header)
    nrdfh.write(header)
    filteredfh.write(header)
    concordancefh.write(header)
    #multifh.write(header)

    samples=vcfobj.getSampleList()

    totalrecords=0

    pattern=';set=(.+)'
    for vrec in vcfobj.yieldVcfRecordwithGenotypes(vcfh):
        if ',' in vrec.getAlt() > 1:
            outstring=vrec.toStringwithGenotypes() + "\n"
            multifh.write(outstring)
            #continue


        
        if 'ReferenceInAll' in vrec.getInfo() and options.includeRef == False:
            continue

        if 'filterIn' in vrec.getInfo():
            outstring=vrec.toStringwithGenotypes() + "\n"
            filteredfh.write(outstring)
            continue
        if 'FilteredInAll' in vrec.getInfo():
            outstring=vrec.toStringwithGenotypes() + "\n"
            filteredfh.write(outstring)
            continue
        vrec_ziptuple=vrec.zipGenotypes(samples)

        #what set are you in?
        field=re.search(pattern, vrec.getInfo()).groups()[0]
        fieldsfh.write(field+"\n")
        totalrecords+=1
        for (compare, eval) in grouper(2,vrec_ziptuple):

           
            (comp_allele1, comp_allele2)=compare[1].getAlleles()
            (eval_allele1, eval_allele2)=eval[1].getAlleles()

            eval_alleletype=typeofGenotype(eval_allele1, eval_allele2)
            comp_alleletype=typeofGenotype(comp_allele1, comp_allele2)
           
            concordancetable[eval_alleletype, comp_alleletype]+=1



            #print records that contirubut the NRS penalty
            if eval_alleletype == 3:
                if comp_alleletype == 1 or comp_alleletype==2:
                    outstring=vrec.toStringwithGenotypes() + "\n"
                    nrsfh.write( outstring)
            if eval_alleletype==0:
                if comp_alleletype == 1 or comp_alleletype == 2:
                    outstring=vrec.toStringwithGenotypes() + "\n"
                    nrsfh.write( outstring )
    
        
            #print records that contribute to NRD penalty
            if eval_alleletype==0:
                if comp_alleletype == 1 or comp_alleletype == 2:
                    outstring=vrec.toStringwithGenotypes() + "\n"
                    nrdfh.write( outstring )
                if comp_alleletype == 0:
                    outstring=vrec.toStringwithGenotypes() + "\n"
                    concordancefh.write( outstring )
            if eval_alleletype == 1:
                if comp_alleletype == 0 or comp_alleletype == 2:
                    outstring=vrec.toStringwithGenotypes() + "\n"
                    nrdfh.write( outstring )
                if comp_alleletype == 1:
                    outstring=vrec.toStringwithGenotypes() + "\n"
                    concordancefh.write( outstring )
            if eval_alleletype == 2:
                if comp_alleletype == 0 or comp_alleletype ==1:
                    outstring=vrec.toStringwithGenotypes() + "\n"
                    nrdfh.write( outstring )
                if comp_alleletype == 2:
                    outstring=vrec.toStringwithGenotypes() + "\n"
                    concordancefh.write( outstring )


    print "total records analyzed: " + str(totalrecords) + "\n"

    print "rows are eval genotypes columns comparison genotypes"
    print "\n"
    print "\t".join(['','AA','AB','BB', './.'  ])
   
    rownames=[0,'AA', 1,'AB', 2,'BB', 3,'./.']
    for (i, gt) in grouper(2,rownames):
        row=concordancetable[i,:].tolist()
        for r in row:
            outstr="\t".join(map(str,r))
            print gt,"\t", outstr

    print "matrix sum: "
    sum=np.sum(concordancetable)
    print str(sum)
    print "\n"

    #now we figure out how many sites were called or not called
    calledtable[0,0]=concordancetable[0:3,0:3].sum()
    calledtable[0,1]=concordancetable[0:3,3].sum()
    calledtable[1,0]=concordancetable[3,0:3].sum()
    calledtable[1,1]=concordancetable[3,3]

    rownames=[ 0,'called', 1,'./.' ]
    print "rows are called eval genotypes columns are called comparison genotypes"
    print "\n"
    print "\t".join(['','called','./.' ])
    
    for (i, gt) in grouper(2,rownames):
        row=calledtable[i,:].tolist()
        for r in row:
            outstr="\t".join(map(str,r))
            print gt,"\t", outstr
    print "matrix sum: "
    sum=np.sum(calledtable)
    print str(sum)
    print "\n"



    if options.matrixonly == False:
        discordance=concordancetable[0,1]+concordancetable[0,2]+concordancetable[1,0]+concordancetable[1,2]+concordancetable[2,0]+concordancetable[2,1]
        total=concordancetable[0,1]+concordancetable[0,2]+concordancetable[1,0]+concordancetable[1,1]+ concordancetable[1,2]+concordancetable[2,0]+concordancetable[2,1] +concordancetable[2,2]
    
        nrd=round( (float(discordance)/float(total)) * 100, 2)
    
        variant_count_evaluation= concordancetable[1,1]+ concordancetable[1,2]+ concordancetable[2,1]+ concordancetable[2,2]
    
        variant_count_comparison= concordancetable[0,1]+concordancetable[0,2]+concordancetable[1,1]+concordancetable[1,2]+concordancetable[2,1]+concordancetable[2,2]+concordancetable[3,1]+concordancetable[3,2]
        nrs=round( float(variant_count_evaluation)/float(variant_count_comparison) * 100 , 2)
    
        print "NRD: ", str(nrd)
        print "NRS ", str(nrs)

   
# <codecell>

if __name__ == "__main__":
    main()

