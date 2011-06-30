from VcfMetaLines import MetaLines
from VcfMetaLines import HeaderLines
from VcfRecord import VcfRecord
from VcfGenotype import VcfGenotype

import string
import sys

class VcfFile(object):
    'a class representing a VCF file'

    def __init__(self, filename):
        self.metaline = MetaLines()
        self.headerline = HeaderLines()
        self.filename=filename


    def parseMetaLines(self, fh):
        """ parse the meta lines that begin with ## in a VCF file """
        for line in fh:
            if '##fileformat' in line.strip():
                self.metaline.setFileFormat(line)
            elif '##INFO' in line.strip():
                self.metaline.parseMetaInfo(line)
            elif '##FILTER' in line.strip():
                self.metaline.parseMetaFilter(line)
            elif '##FORMAT' in line.strip():
                self.metaline.parseMetaFormat(line)
            elif '#CHROM' in line.strip():
                break
        

    def hasFormatHeader(self,headerline):
        """ check if header line has FORMAT column """
        self.headerline.add_format_column(headerline)


    def parseHeaderLine(self,fh):
        for line in fh:
            if '#CHROM' not in line: continue
            
            self.headerline.add_format_column(line.strip() )
            self.headerline.append_samplelist( line.strip() )
            
            break
            

    

    def printHeaderLine(self):
        print self.headerline.toString()



    def getMetaInfoNames(self):
        """  return the IDs of the INFO metalines"""
        return self.metaline.getMetaInfoNames()

    def getMetaInfoDescription(self):
        return self.metaline.getMetaInfoDescription()

    def getMetaFilterNames(self):
        """ return the IDs of the FILTER metalines """
        return self.metaline.getMetaFilterNames()
    
    def getMetaFilterDescription(self):
        """ return a list of tuples with (id, description) of FILTER metalines """
        return self.metaline.getMetaFilterDescription()

    def getMetaFormatNames(self):
        """  return the Ids of the FORMAT metalines"""
        return self.metaline.getMetaFormatNames()

    def getMetaFormatDescription(self):
        """ return a list of tuples (id,descripton) for FORMAt metalines """
        return self.metaline.getMetaFormatDescription()

    def yieldMetaInfoLines(self):
        """ yield string representation of  ##INFO metalines """
        return self.metaline.yieldPrintMetaInfoLines()

    def yieldMetaFormatLines(self):
        """ yield string reprsentation ##FORMAT metaline object """
        return self.metaline.yieldPrintMetaFormatLines()

    def yieldMetaFilterLines(self):
        """ yield string represenation of ##FILTER metaline object """
        return self.metaline.yieldPrintMetaFilterLines()


    def printMetaInfoLines(self):
        """print the meta ##INFO lines for a VCF """
        for str in self.yieldMetaInfoLines():
            print str

    def printMetaFormatLines(self):
        for str in self.yieldMetaFormatLines():
            print str
    
    def printMetaFilterLines(self):
        for str in self.yieldMetaFilterLines():
            print str

    def printMetaLines(self):
        self.printMetaInfoLines()
        self.printMetaFormatLines()
        self.printMetaFilterLines()

    def yieldDatalines(self, fh):
        """ yield a dataline from a VCF file """
        for line in fh:
            yield line.strip()

    def yieldVcfRecord(self,fh):
        """ yield VcfRecord object from reading a dataline in a VCF file """
        for line in fh:
            if '#' in line: continue
            fields=line.strip().split('\t')
            (chrom,pos,id,ref,alt,qual,filter,info)=fields[0:8]
            yield VcfRecord(chrom,pos,id,ref,alt,qual,filter,info)

    def yieldVcfRecordwithGenotypes(self,fh):
        """ yield a VcfRecord with its genotypes list populated """
        for line in fh:
            fields=line.strip().split('\t')
            (chrom,pos,id,ref,alt,qual,filter,info)=fields[0:8]
            vrec=VcfRecord(chrom,pos,id,ref,alt,qual,filter,info)
    
            formatstring=fields[8]
            (genotypestrings)=fields[9::]
            for gstring in genotypestrings:
                vrec.addGenotype( VcfGenotype(formatstring,gstring))
            yield vrec

    def printDataLineWithGenotypes(self,vcfh):
        for vrec in self.yieldVcfRecordwithGenotypes(vcfh):
            print vrec.toStringwithGenotypes()