import sys
import re
class VcfGenotype(object):
    """ represents a VCF genotype  """

    def __init__(self, formatstring, gstring ):
        """ initialize a VcfGenotype """
        """ formatstring is from the format column gstring is the genotypestring """
        self.gstring=gstring
        self.formatstring=formatstring
        self.isPhased=0
        self.allele1=''
        self.allele2=''

        #keys are format key value is from gstring
        self.gdict={}
        formatids=self.formatstring.split(':')
        gstringvals=self.gstring.split(':')
        self.parseAlleles(gstringvals[0])

        if gstring != '.':
            if len(formatids) != len(gstringvals):
                sys.stderr.write("error mismatch btwn format and genotype string in VcfGenotype init!\n")
                print len(formatids)
                print len(gstringvals)
                exit(1)

            for i in range(0, len(formatids) ):
                self.gdict[ formatids[i] ] = gstringvals[i]
        else:
            for i in range(0, len(formatids) ):
                self.gdict[ formatids[i] ] = '.'

    def getGenotypeFormatFields(self):
        """ return list of ids for genotype format string """
        return self.gdict.keys()


    def setAlleles(self, allele1, allele2):
        """ set the alleles for the VcfGenotype object """
        self.allele1=allele1
        self.allele2=allele2

    def setIsPhased(self):
        """ set phased flag on the gentoype to true """
        self.isPhased=1

    def parseAlleles(self, allelefield):
        """ set the allele1 and allele2 given a genotype (GT) field """
        delimiter=''
        if '|' in allelefield:
            setIsPhased()
            delimiter='|'
        elif '/' in allelefield:
            delimiter='/'
        elif '.' == allelefield:
            self.allele1='.'
            self.allele2='.'
            return
        else:
            sys.stderr.write("un-recognized genotype delimiter: " + allelefield + "\n")
            exit(1)
        (allele1,allele2) = allelefield.split(delimiter)

        self.allele1=allele1
        self.allele2=allele2

    def getFormatString(self):
        return self.formatstring

    def getAlleles(self):
        """ return tuple with (allele1, allele2) """
        return (self.allele1, self.allele2)

    def isCalled(self):
        """ return True if genotype was called and both alleles are not '.' """
        if self.allele1 != '.' and self.allele2 != '.':
            return True
        return False

    def isHet(self):
        """ return True if gentoype is het and both alleles are called"""
        
        if self.allele1 != '.' and self.allele2 != '.':
            if self.allele1 != self.allele2:
                return True
        return False
    def checkFormatIds(self, formatlist):
        """ check the ids in the FORMAT field to make sure they are contained in the list of FORMAT ids (the format list) """
        formatfields=self.formatstring.split(':')
        for elem in formatfields:
            if elem not in self.gdict.keys():
                sys.stderr.write(elem + " is not in FORMAT column!\n ")
                exit(1)



    def toString(self):
        return self.gstring
