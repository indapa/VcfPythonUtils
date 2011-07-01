import re
import sys
class VcfRecord(object):
    """ reprsents a VCF dataline """

    def __init__(self, chrom='.', pos='.', id='.', ref='.', alt='.', qual='.', filter='.', info='.'):
        """ iniialize a VCF record """
        self.chrom=chrom
        self.pos=pos
        self.id=id
        self.ref=ref
        self.alt=alt
        self.qual=qual
        self.filter=filter
        self.info=info
       
        self.genotypes=[] # list of VcfGenotype object as elements

    def setChrom(self,chrom):
        self.chrom=chrom

    def setPos(self,pos):
        self.pos=pos

    def setId(self,id):
        self.id=id

    def setRef(self,ref):
        self.ref=ref

    def setAlt(self,alt):
        self.alt=alt

    def setQual(self,qual):
        self.qual=qual

    def setFilter(self,filter):
        self.filter=filter

    def setInfo(self,info):
        self.info=info



    def getChrom(self):
        return self.chrom

    def getPos(self):
        return self.pos

    def getId(self):
        return self.id

    def getRef(self):
        return self.ref

    def getAlt(self):
        return self.alt

    def getQual(self):
        return self.qual

    def getFilter(self):
        return self.filter

    def getFilter(self):
        return self.info

    def checkInfoIds(self, infolist):
        """ check the ids in the INFO field to make sure they are contained in the list of INFO ids ( the infolist)  """
        pattern= '(.*)='
        infofields=self.info.split(';')
        for elem in infofields:
            if '=' not in elem:
                id = elem
            else:
                id =re.search(pattern, elem).groups()[0]
            if id == None:
                sys.stderr.write("error in parsing INFO column in VcfRecord!\n")
                exit(1)
            else:
                if id not in infolist:
                    sys.stderr.write(id + " not in ##INFO header!\n")
                    exit(1)


    def addGenotype(self,genotypeobj):
        """ append a VcfGenotype object to genotype list """
        self.genotypes.append(genotypeobj)

    def toString(self):
        outstring="\t".join([self.chrom,self.pos,self.id,self.ref,self.alt,self.qual,self.filter,self.info])
        return outstring

    def toStringwithGenotypes(self):
        outstring=self.toString()
        formatstring=self.genotypes[0].getFormatString()

        genotypestrings=[]
        for g in self.genotypes:
            genotypestrings.append( g.toString() )

        genostring="\t".join(genotypestrings)
        return outstring + "\t" + formatstring + "\t" + genostring

