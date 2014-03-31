'''
Created on May 28, 2013

@author: indapa
'''
from VcfFile import *
import numpy as np
from common import grouper, melt_lol

class VcfSampleEval(object):
    """ a class to represent a genotype comparison 
        We can do a sample by sample comparison, rather than
        lump all samples together when evalating a callset like in variantEvalGenotypeConcordance.py"""
    def __init__(self, compareName, evalName, basename, logfiles=False):
        
        self.comparename=compareName
        self.evalname=evalName
        
        
        self.concordancetable= np.matrix( [ [ 0,0,0,0 ], [ 0,0,0,0 ], [ 0,0,0,0 ], [ 0,0,0,0 ] ] )
        self.calledtable = np.matrix ( [ [0 ,0] , [0,0] ] )
        
        self.outputfile=".".join([basename, evalName, compareName,'variantEval','txt'])
        self.nrslog=".".join([basename, evalName, compareName,'nrs','log'])
        self.nrdlog=".".join([basename, evalName, compareName,'nrd','log'])
        self.concordancelog=".".join([basename,evalName, compareName, 'concordance','log'])
        self.genotypematrix=".".join([basename,evalName, compareName, 'genotype.matrix', 'csv'])
        
        if logfiles == True:
            self.nrsfh=open(self.nrslog, 'w')
            self.nrdfh=open(self.nrdlog, 'w')
            self.concordancefh=open(self.concordancelog, 'w')
            self.outputfh=open(self.outputfile, 'w')
        
        self.genotypematrixfh=open(self.genotypematrix, 'w')
        
        
        
    def write_genotype_matrix(self):
        """ melt the genotypematrix into a CSV of numbers  """
        outstring=",".join( map(str,melt_lol(self.concordancetable.tolist())) )
        self.genotypematrixfh.write(outstring+"\n")
        
    def incrementcellcount (self, eval_alleletype, comp_alleletype):
        """ the parameters are results of calling common.typeofGenotype
         which results in returning 0 (homoz_ref) 1 het 2 homoz_alt or 3 nocall
         then just increment the appropriate row,column in the concordance matrix
         """
        self.concordancetable[eval_alleletype, comp_alleletype]+=1
        
    
    def writeHeaders(self, headerline):
        """ write VCF headers to the log files  """
        self.nrdfh.write( headerline )
        self.nrsfh.write( headerline )
        self.concordancefh.write( headerline )

    def writeNrd(self,outstring):
        """ write the vcfrecord dataline outstirng to the nrd log """
        self.nrdfh.write( outstring )
         
    def writeNrs(self,outstring):
        """ write the vcfrecord dataline outstirng to the nrs log """
        self.nrsfh.write( outstring )
    
    
    def writeConcordance(self,outstring):
        """ write the vcfrecord dataline outstirng to the concordance log """
        self.concordancefh.write( outstring )
        
    
    def returnNRS_NRD(self):
        
        discordance=self.concordancetable[0,1]+self.concordancetable[0,2]+self.concordancetable[1,0]+self.concordancetable[1,2]+self.concordancetable[2,0]+self.concordancetable[2,1]
        total=self.concordancetable[0,1]+self.concordancetable[0,2]+self.concordancetable[1,0]+self.concordancetable[1,1]+ self.concordancetable[1,2]+self.concordancetable[2,0]+self.concordancetable[2,1] +self.concordancetable[2,2]
    
        nrd=round( (float(discordance)/float(total)) * 100, 2)
    
        variant_count_evaluation= self.concordancetable[1,1]+ self.concordancetable[1,2]+ self.concordancetable[2,1]+ self.concordancetable[2,2]
    
        variant_count_comparison= self.concordancetable[0,1]+self.concordancetable[0,2]+self.concordancetable[1,1]+self.concordancetable[1,2]+self.concordancetable[2,1]+self.concordancetable[2,2]+self.concordancetable[3,1]+self.concordancetable[3,2]
        nrs=round( float(variant_count_evaluation)/float(variant_count_comparison) * 100 , 2)
        
        return (nrs, nrd)
    
    def writeEvalOutput(self):
        
        self.outputfh.write( "rows are eval genotypes columns comparison genotypes\n")
    
        self.outputfh.write("\t".join(['','AA','AB','BB', './.'  ])  +"\n")
   
        rownames=[0,'AA', 1,'AB', 2,'BB', 3,'./.']
        for (i, gt) in grouper(2,rownames):
            row=self.concordancetable[i,:].tolist()
            for r in row:
                outstr="\t".join(map(str,r))
                self.outputfh.write( gt +"\t"+outstr+"\n")

        self.outputfh.write( "matrix sum: \n")
        summy=np.sum(self.concordancetable)
        self.outputfh.write( str(summy) +"\n")
        
        
        
        #now we figure out how many sites were called or not called
        self.calledtable[0,0]=self.concordancetable[0:3,0:3].sum()
        self.calledtable[0,1]=self.concordancetable[0:3,3].sum()
        self.calledtable[1,0]=self.concordancetable[3,0:3].sum()
        self.calledtable[1,1]=self.concordancetable[3,3]
        self.outputfh.write("\n")
        rownames=[ 0,'called', 1,'./.' ]
        self.outputfh.write( "rows are eval genotypes columns comparison genotypes\n")
    
        self.outputfh.write(  "\t".join(['','called','./.' ]) +"\n" )
    
        for (i, gt) in grouper(2,rownames):
            row=self.calledtable[i,:].tolist()
            for r in row:
                outstr="\t".join(map(str,r))
                self.outputfh.write( gt +"\t"+outstr+"\n")
        self.outputfh.write( "matrix sum: \n")
        summy=np.sum(self.calledtable)
        self.outputfh.write( str(summy) +"\n")
   
        self.outputfh.write("\n")
        
        
        discordance=self.concordancetable[0,1]+self.concordancetable[0,2]+self.concordancetable[1,0]+self.concordancetable[1,2]+self.concordancetable[2,0]+self.concordancetable[2,1]
        total=self.concordancetable[0,1]+self.concordancetable[0,2]+self.concordancetable[1,0]+self.concordancetable[1,1]+ self.concordancetable[1,2]+self.concordancetable[2,0]+self.concordancetable[2,1] +self.concordancetable[2,2]
    
        nrd=round( (float(discordance)/float(total)) * 100, 2)
    
        variant_count_evaluation= self.concordancetable[1,1]+ self.concordancetable[1,2]+ self.concordancetable[2,1]+ self.concordancetable[2,2]
    
        variant_count_comparison= self.concordancetable[0,1]+self.concordancetable[0,2]+self.concordancetable[1,1]+self.concordancetable[1,2]+self.concordancetable[2,1]+self.concordancetable[2,2]+self.concordancetable[3,1]+self.concordancetable[3,2]
        nrs=round( float(variant_count_evaluation)/float(variant_count_comparison) * 100 , 2)
    
        self.outputfh.write( "NRD: " + str(nrd) +" \n")
        self.outputfh.write( "NRS " + str(nrs) +" \n")
        
        outstring=",".join( map(str,melt_lol(self.concordancetable.tolist())) )
        self.genotypematrixfh.write(outstring+"\n")
        
         
        
    def __str__(self):
        rownames=[0,'AA', 1,'AB', 2,'BB', 3,'./.']
        outstring="\t".join(['','AA','AB','BB', './.']) +"\n" 
        
        for (i, gt) in grouper(2,rownames):
            row=self.concordancetable[i,:].tolist()
            for r in row:
                outstr="\t".join(map(str,r))
                outstring+=( gt +"\t"+outstr+"\n")
        outstring+="eval: "+ self.evalname
        outstring+=" compare: "+ self.comparename +"\n"
        return outstring
    
        
        
        
        
    
    
    