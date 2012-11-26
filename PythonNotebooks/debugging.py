# -*- coding: utf-8 -*-
# <nbformat>2</nbformat>

# <codecell>

from VcfFile import *
vcfh=open('filter.vcf','r')
vcfobj=VcfFile('filter.vcf')
vcfobj.parseMetaAndHeaderLines(vcfh)
header=vcfobj.returnHeader() +"\n"
samplelist=vcfobj.getSampleList()
affecteds=[] # list of affected samples
unaffecteds=[] # list of unaffected samples
affectedfh=open('affecteds.txt', 'r')
for line in affectedfh:
    affecteds.append(line.strip() )
    
unaffectedfh=open ('unaffecteds.txt', 'r')
for line in unaffectedfh:
    unaffecteds.append ( line.strip() )

if len( list( set(unaffecteds).intersection( set(affecteds) ) )  ) != 0:
    sys.stderr.write("check list of affected and unaffecteds for overlapping samples!\n")
    exit(1)
#print samplelist
filter='PASS'
model='recessive'
for vrec in vcfobj.yieldVcfRecordwithGenotypes(vcfh ):

    affected_genotypes=[] ##list of tuples (sample, VcfGenotype object) with samples that are affected
    unaffected_genotypes=[]##list of tuples (sample, VcfGenotype object) with samples that are unaffected
    
    if vrec.getFilter() != filter and filter != None : continue
    #print vrec.toString()
    genotypes = vrec.getGenotypesAlleles()
    genotype_tuple= vrec.zipGenotypes(samplelist) # get a list of tuples [ (sample, VcfGenotype object) ... ]
    for (sample, genotype) in genotype_tuple: #iterate thru and see if they are in affected or unaffected list
        if model == 'dominant':
            if sample in affecteds:  # if so ...
                affected_genotypes.append( ( sample, genotype.toString(),  genotype.isSegregating() )  ) # are they segregating for a non-ref allele?
            if sample in unaffecteds:
                unaffected_genotypes.append( (sample,  genotype.toString(),  genotype.isSegregating() ) ) # are they segregating for a non-ref allele?
        elif model == 'recessive':
            if sample in affecteds:
                affected_genotypes.append( ( sample, genotype.toString(),  genotype.isNonRefHomz() )  ) # are they segregating for a non-ref homoz?
            if sample in unaffecteds:
                unaffected_genotypes.append( (sample,  genotype.toString(),  genotype.isNonRefHomz() ) ) # are they segregating for a non-ref non-refhomoz?
        else:
            sys.stderr.write(model + " not supported for genotype discrete filtering ...\n")
#    print "\n"
#    print affected_genotypes
#    print "\n"
#    print unaffected_genotypes

    shared_affected_segregating=[]

    shared_unaffected_segregating=[]
    if model == 'dominant':
        #print 'under dominant inhertiance ....\n'
        """ under dominant model, all affecteds should be segrgating for non-ref allele
            and all UN-affecteds should *NOT* be segregating for non-ref allele  """
 
        
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
    elif model == 'recessive':

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
        sys.stderr.write(model + " not supported for genotype discrete filtering ...\n")
# <codecell>


