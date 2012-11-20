class Ped(object):
    """ represents a ped file http://pngu.mgh.harvard.edu/~purcell/plink/data.shtml#ped """

    def __init__(self,famid='.', indv='.', paternal='.', maternal='.', sex='.', phenotype='.'):
        """ attributes of a Ped object """
        self.famid=famid
        self.individ=indv
        self.pid=paternal
        self.mid=maternal
        self.sex=sex
        self.pheno=phenotype


    def setfamid(self,famid):
        self.famid=famid

    def setindvid(self,indv):
        self.individ=indv

    def setmid(self,mid):
         self.mid=mid

    def setpid(self,pid):
         self.pid=pid

    def setsex(self,sex):
         self.sex=sex

    def setpheno(self,pheno):
         self.pheno=pheno



    def getfamid(self):
        return self.famid

    def getid(self):
        return self.individ

    def getmid(self):
        return self.mid

    def getpid(self):
        return self.pid

    def getsex(self):
        return self.sex

    def getpheno(self):
        return self.pheno

    def toString(self):
        return "\t".join( [ self.famid, self.individ, self.pid, self.mid, self.sex, self.pheno] )
    
