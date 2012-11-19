#!/usr/bin/env python

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

#for genotype_space in doCrossAutosomal( ('.','.'), ('.', '.') ):
for genotype_space in doCrossAutosomal( ('0','0'), ('.', '1') ):
    print genotype_space

alleles=('0','1')
maternal=('.', '1')
paternal=('0','1')

