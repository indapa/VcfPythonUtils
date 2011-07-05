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

    def getInfo(self):
        return self.info


    def check_genotypeFormat(self, formatlist):
        """ check the ids in the FORMAT column are contained in the list of FORMAT ids (formatlist) """
        if len(self.genotypes) == 0:
            sys.stderr.write("VCF file contains no genotype columns\n")
            return
        else:
            self.genotypes[0].checkFormatIds(formatlist)

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

    def appendInfoString(self, infostring ):
        """ apppend to INFO string for a VCF record """
        self.info+=";"+infostring

    def addGenotype(self,genotypeobj):
        """ append a VcfGenotype object to genotype list """
        self.genotypes.append(genotypeobj)


    def getGenotypesAlleles(self):
        """ return list of tuples with (allele1, allele2) for each VcfGenotype object in genotypes list  """
        genotypeAlleles=[]
        for vcfgenobj in self.genotypes:
            genotypeAlleles.append( vcfgenobj.getAlleles() )
       
        return genotypeAlleles

    def allHets(self):
        """ return true if all genotypes for the VcfREcord are hets """
        for gobj in self.genotypes:
            if gobj.isHet() == False:
                return False
        return True

    def siteCallrate(self):
        """ calculate the site callrate: #called_genotypes/#genotypes """
        called=0
        for gobj in self.genotypes:
            if gobj.isCalled() == True:
                called+=1
        callrate=float(called)/float(len(self.genotypes))
        return callrate

    def sampleCallrate(self,samplist,samplecallsdict):
        called=[]
        ziplist= zip(samplist, self.genotypes)
        for (sample, gobj) in ziplist:
            if  gobj.isCalled() == True:
                samplecallsdict[sample]+=1
            

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




