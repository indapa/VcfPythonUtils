import sys
import re
import itertools
from collections import OrderedDict
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
        self.gdict=OrderedDict()
        formatids=self.formatstring.split(':')
        gstringvals=self.gstring.split(':')
        
        self.parseAlleles(gstringvals[0]) # the first format is always the GT (genotype)
       
        zipiter=itertools.izip_longest(formatids,gstringvals,fillvalue='.')

        for (format,gstringval) in zipiter:
            self.gdict[format]=gstringval

       

    def getGenotypeFormatFields(self):
        """ return list of ids for genotype format string """
        return self.gdict.keys()


    def getFormatVal(self, key):
        """ return a value from the format string  """
        if key in self.gdict:
            return self.gdict[key]
        else:
            return "."
        
    def setFormatVal(self,key,value):
        """ set a value of a format field """
        if key in self.gdict:
            self.gdict[key]=value
            """ then reset the gstring"""
            self.gstring=":".join( [ self.gdict[k] for k in self.gdict.keys() ] )
        else:
            sys.stderr.write(key + " not in format of genotype!\n")
            
        
        
    def addFormatVal(self, key, value):
        """ add a new format key/value to the genotype format string;
            we then implicitly update the gstring and formatstring  """
       
        self.gdict[key]=value
        self.gstring=":".join( [ self.gdict[k] for k in self.gdict.keys() ] )
        self.formatstring=":".join(self.gdict.keys())
        
    def addFormat(self,formatstring):
        """ add a format field to the genotype """ 
        self.formatstring+=":"+formatstring
         
        

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
            self.setIsPhased()
            delimiter='|'
        elif '/' in allelefield:
            delimiter='/'
        elif '.' == allelefield:
            self.allele1='.'
            self.allele2='.'
            return
        else:
            sys.stderr.write("un-recognized genotype delimiter: " + allelefield + "\n")
            self.allele1='.'
            self.allele2='.'
            return
            
        (allele1,allele2) = allelefield.split(delimiter)

        self.allele1=allele1
        self.allele2=allele2

    def setFormatAndGString(self, fstring, gstring):
        #print gstring
        self.formatstring=fstring
        self.gstring=gstring
        
        self.gdict=OrderedDict()
        formatids=self.formatstring.split(':')
        gstringvals=self.gstring.split(':')
        
        self.parseAlleles(gstringvals[0]) # the first format is always the GT (genotype)
       
        zipiter=itertools.izip_longest(formatids,gstringvals,fillvalue='.')

        for (format,gstringval) in zipiter:
            self.gdict[format]=gstringval
        #print self.gstring
        
    
    def setFormatString(self,fstring):
        self.formatstring =fstring
        
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


    def isSegregating(self):
        """ return true if genotype has at least one non-ref allele so 1/1 and 1/0 would return true but 0/0 returns false """
        if self.isCalled() == False:
            return False
        if self.allele1 != '0' or self.allele2 != '0':
            return True
        return False

    def isNonRefHomz(self):
        """ return true if genotype is non-ref homozygote """
        if self.isCalled() == False:
            return False
        if self.allele1 !='0' and self.allele2 != '0':
            return True
        return False



    def isHet(self):
        """ return True if gentoype is het and both alleles are called"""
        
        if self.allele1 != '.' and self.allele2 != '.':
            if self.allele1 != self.allele2:
                return True
        return False
    
    def isHomoz(self):
        if self.allele1 != "." and self.allele2 != ".":
            if self.allele1 == '0' and self.allele2 == "0":
                return True
        return False
    
    def checkFormatIds(self, formatlist):
        """ check the ids in the FORMAT field to make sure they are contained in the list of FORMAT ids (the format list) """
        formatfields=self.formatstring.split(':')
        for elem in formatfields:
            if elem not in self.gdict.keys():
                sys.stderr.write(elem + " is not in FORMAT column!\n ")
                exit(1)
    def containsFormatKey(self, key):
        if key in self.gdict.keys():
            return True
        return False
    
    def toString(self):
        return self.gstring

    def __str__(self):
        return self.gstring

