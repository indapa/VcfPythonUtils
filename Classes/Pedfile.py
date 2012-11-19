from VcfPed import *

class Pedfile(object):
    """ a Pedfile object has a list of VcfPed objects """

    def __init__(self,filename):
        self.filename=filename
        self.vcfpedlist=[]

    def parsePedfile(self, fh):
        """ given a filehandle to a *.ped file read its contents and populate the list vcfpedlist with VcfPed objects """
        for line in fh:
            fields=line.strip().split('\t')
            (famid, indv, pid, mid, sex, pheno)=fields[0:6]
            self.vcfpedlist.append( Ped(famid, indv, pid, mid, sex, pheno) )

    def toString(self):
        for obj in self.vcfpedlist:
            print obj.toString()

    def returnFounders(self):
        """ return the founders in a ped file (those with unknown paternal and maternids """
        founders=[]

        for vcfpedobj in self.vcfpedlist:
            if vcfpedobj.getpid() == "0":
                founders.append(vcfpedobj)

        return founders

    def returnFounderIds(self):
        """ return the indiv ids of the founders in the ped file"""
        founderids=[]
        for vcfpedobj in self.vcfpedlist:
            if vcfpedobj.getpid() == "0":
                founderids.append( vcfpedobj.getid() )
        return founderids

    def returnNonFounderIds(self):
        """ return the indiv ids of the founders in the ped file"""
        nonfounderids=[]
        for vcfpedobj in self.vcfpedlist:
            if vcfpedobj.getpid() != "0":
                nonfounderids.append( vcfpedobj.getid() )
        return nonfounderids

    def returnNonFounders(self):
        """ return the founders in a ped file (those with unknown paternal and maternids """
        nonfounders=[]

        for vcfpedobj in self.vcfpedlist:
            if vcfpedobj.getpid() != "0":
                nonfounders.append(vcfpedobj)

        return nonfounders

    def returnIndivids(self):
        """ return a list of indvi ids from the ped file """
        #samplelist=[]
        return [ ped.getid() for ped in self.vcfpedlist]
