#!/usr/bin/env python
import sys
import gzip
from optparse import OptionParser

from VcfFile import *
from VcfPed import Ped

""" Given a ped file with affected and unaffected status in the 6th column
    and a given inheritance model ( dominant|recessive)
    filter genotypes in a VCF to produce a sites as
    possible candidates for causative mutation for
    a Mendelian trait. For more on the pedfile fomrat see this: http://pngu.mgh.harvard.edu/~purcell/plink/data.shtml#ped"""


def main():
    usage = "usage: %prog [options] file.vcf.gz"
    parser = OptionParser(usage)
    parser.add_option("--model", type="string", dest="model", default = "dominant", help=" inheritance model [dominant|recessive], default is dominant ")
    parser.add_option("--ped", type="string", dest="pedfile", default=None, help="ped file of samples with phenotype (disease) status")
    parser.add_option("--filter", type="string", dest="filter", help="analyze only those  records matching filter (default is PASS)", default='PASS')

    (options, args)=parser.parse_args()
    if options.pedfile==None:
        sys.stderr.write("please provide a value to --ped parameter!\n")
        exit(1)


    affecteds=[] # list of affected samples
    unaffecteds=[] # list of unaffected samples
    
    pedobjects=[] #list of pedobjects, represents lines in a pedfile
    pedfh=open(options.pedfile, 'r')
    for line in pedfh:
        fields=line.strip().split('\t')
        (fid,iid,pid,mid,sex,phenotype)=fields[0:6]
        phenotype=int(phenotype)
        pedobjects.append( Ped(fid,iid,pid,mid,sex,phenotype) )

    #the phenotype status is set to 2 if the sample is affected: http://pngu.mgh.harvard.edu/~purcell/plink/data.shtml#ped
    affecteds=[ pedobj.getid() for pedobj in pedobjects if pedobj.getpheno() == 2  ]
    unaffecteds=[ pedobj.getid() for pedobj in pedobjects if pedobj.getpheno() == 1  ]



    

    #check if any overlapping samples between unaffected and affected
    if len( list( set(unaffecteds).intersection( set(affecteds) ) )  ) != 0:
        sys.stderr.write("check list of affected and unaffecteds for overlapping samples!\n")
        exit(1)

    #    sys.stderr.write("check list of affected and unaffected for overlapping samples!\n")
    #    exit(1)


    vcfilename=args[0]
    vcfh=gzip.open(vcfilename,'r')

    #instantiate a VcfFile object
    vcfobj=VcfFile(vcfilename)
    vcfobj.parseMetaAndHeaderLines(vcfh)
    header=vcfobj.returnHeader()
    samplelist=vcfobj.getSampleList()

    print header

    for vrec in vcfobj.yieldVcfRecordwithGenotypes(vcfh ):
        
        affected_genotypes=[] #list of tuples (sample, VcfGenotype object) with samples that are affected
        unaffected_genotypes=[] # list of tuples (sample, VcfGenotype object) with samples that are unaffected

        if vrec.getFilter() != options.filter and options.filter != None : continue
        
        genotype_tuple= vrec.zipGenotypes(samplelist) # get a list of tuples [ (sample, VcfGenotype object) ... ]
        for (sample, genotype) in genotype_tuple: #iterate thru and see if they are in affected or unaffected list
            if options.model == 'dominant':
                if sample in affecteds:  # if so ...
                    affected_genotypes.append( ( sample, genotype.toString(),  genotype.isSegregating() )  ) # are they segregating for a non-ref allele?
                if sample in unaffecteds:
                    unaffected_genotypes.append( (sample,  genotype.toString(),  genotype.isSegregating() ) ) # are they segregating for a non-ref allele?
            elif options.model == 'recessive':
                if sample in affecteds:
                    affected_genotypes.append( ( sample, genotype.toString(),  genotype.isNonRefHomz() )  ) # are they segregating for a non-ref homoz?
                if sample in unaffecteds:
                    unaffected_genotypes.append( (sample,  genotype.toString(),  genotype.isNonRefHomz() ) ) # are they segregating for a non-ref non-refhomoz?
            else:
                sys.stderr.write(options.model + " not supported for genotype discrete filtering ...\n")


        if options.model == 'dominant':
        #under dominant model, all affecteds should be
        #segrgating for non-ref allele and all UN-affecteds should *NOT* be segregating for non-ref allele
            
            #how many affected individuals are segregating for non-ref allele?
            count_segregating_affected = [ tpl[2] == True for tpl in affected_genotypes ].count(True)

            #how many UN-affected individuals are *NOT*  segregating for non-ref allele?
            count_segregating_unaffected =  [ tpl[2] == False for tpl in unaffected_genotypes ].count(True)

            #now if all affects are segregating for the site
            # and all the un-affecteds are *not* segregating for the site
            # it is a candidate
            if count_segregating_affected == len(affecteds):
                if  count_segregating_unaffected == len(unaffecteds):
                    print vrec.toStringwithGenotypes()

        elif options.model == 'recessive':
            #how many affected individuals are segregating for non-ref allele?
            #http://stackoverflow.com/a/5684324/1735942
            count_homoz_nonref_affected = [ tpl[2] == True for tpl in affected_genotypes ].count(True)

            #how many UN-affected individuals are *NOT*  segregating for non-ref allele?
            count_homoz_ref_unaffected =  [ tpl[2] == False for tpl in unaffected_genotypes ].count(True)



            #now if all affects are homoz nonref for the site
            # and all the un-affecteds are homoz ref for the site
            # it is a candidate
            if count_homoz_nonref_affected == len(affecteds):
                if  count_homoz_ref_unaffected  == len(unaffecteds):
                    print vrec.toStringwithGenotypes()
        else:
            sys.stderr.write(options.model + " not supported for genotype discrete filtering ...\n")



if __name__ == "__main__":
    main()
