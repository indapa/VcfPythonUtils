import re
import sys

class InfoLine(object):
    """ represents an ##INFO line  """
    """ ID=ID,Number=number,Type=type,Description=description """
    def __init__(self, id='', number='', type='', description=''):
        self.id=id
        self.number=number
        self.type=type
        self.description=description

    def setId(self,id):
        self.id=id
    def getId(self):
        return self.id

    def setNumber(self, number):
        self.number=number
    def getNumber(self):
        return self.number

    def setType(self,type):
        self.type=type
    def getType(self):
        return self.type

    def setDescription(self,description):
        self.description=description
    def getDescription(self):
        return self.description

    def toString(self):
        """ return string representation of INFO line object """
        infostring=",".join(["ID="+self.id, "Number="+self.number, "Type="+self.type, "Description="+"\""+ self.description + "\""] )
        hashstring="##INFO=<"
        outstring=hashstring+infostring+">"
        return outstring

class FormatLine(object):
    """ represents a ##FORMAT line'"""
    """ ##FORMAT=<ID=ID,Number=number,Type=type,Description=description> """
    """ need to figure out how python inheritance works ..."""
    def __init__(self, id='', number='', type='', description=''):
        self.id=id
        self.number=number
        self.type=type
        self.description=description

    def setId(self,id):
        self.id=id
    def getId(self):
        return self.id

    def setNumber(self, number):
        self.number=number
    def getNumber(self):
        return self.number

    def setType(self,type):
        self.type=type
    def getType(self):
        return self.type

    def setDescription(self,description):
        self.description=description
    def getDescription(self):
        return self.description

    def toString(self):
        """ return string representation of FORMAT line object """

        formatstring=",".join(["ID="+self.id, "Number="+self.number, "Type="+self.type, "Description="+"\""+ self.description + "\""])
        hashstring="##FORMAT=<"
        outstring=hashstring+formatstring+">"
        return outstring

class FilterLine(object):
    """ reprsents a ##FILTER line """
    """ ##FILTER=<ID=ID,Description=description> """
    def __init__(self, id='', description=''):
        self.id=id
        self.description=description

    def setId(self,id):
        self.id=id
    def getId(self):
        return self.id

    def setDescription(self,description):
        self.description=description
    def getDescription(self):
        return self.description

    def toString(self):
        """ return a string representaiton of FILTER line object """
        filterstring=",".join ( ["ID="+self.id, "Description="+ "\""+self.description + "\""] )
        hashstring="##FILTER=<"
        outstring=hashstring+filterstring+">"
        return outstring

class MetaLines(object):

    def __init__(self):
        """ class to hold meta line info from VCF files """
        self.metaFormatDict={}
        """    dict for ##FORMAT lines """

        self.metaFilterDict={}
        """ dictionary for ##FILTER lines """

        self.metaInfoDict={}
        """ dictory for ##INFO lines  """

        self.fileFormat=''


    def getMetaInfoNames(self):
        """ return the INFO IDs in the metaInfoDict """
        return self.metaInfoDict.keys()


    def addMetaInfo(self, id, type, number, description):
        """ add InfoLIne object to metaInfoDict dictionary """
        if id in self.metaInfoDict.keys():
            sys.stderr.write(id + " ##INFO id is already being used!\n")
            exit(1)
            return
        else:
            newinfobject=InfoLine()
            newinfobject.setId(id)
            newinfobject.setType(type)
            newinfobject.setNumber(str(number))
            newinfobject.setDescription(description)
            self.metaInfoDict[ newinfobject.getId() ] = newinfobject
        return

    def addMetaFilter(self, id, description):
        """ add FilterLine object to metaFilterDict dictionary """
        if id in self.metaFilterDict.keys():
            sys.stderr.write(id + "##FILTER id already being used!\n")
            exit(1)
            return
        else:
            newfilterobject=FilterLine(id, description)
            self.metaFilterDict [ newfilterobject.getId() ] = newfilterobject
            
    
    def addMetaFormat(self, formatObj):
        """ add a ##FORMAT header  """
        if formatObj.getId() in self.metaFormatDict.keys():
            sys.stderr.write(id + "##FORMAT id already being used!\n")
            return
        else:
            self.metaFormatDict [ formatObj.getId() ]  = formatObj
            
    



    def getMetaFilterNames(self):
        """ return the FILTER IDs in the metaFilterDict """
        return self.metaFilterDict.keys()

    def getMetaFormatNames(self):
        """ return the FORMAT IDs in metaFormatDict """
        return self.metaFormatDict.keys()


    def getMetaInfoDescription(self):
        """ return a list of tuples with (id,description) of INFO metalines """
        description_strings=[]
        for k in self.metaInfoDict.keys():
            description_strings.append( (k, self.metaInfoDict[k].getDescription() ) )
        return description_strings


    def getMetaFormatDescription(self):
        """ return a list of tuples with (id,description) of FORMAT  metalines """
        description_strings=[]
        for k in self.metaFormatDict.keys():
            description_strings.append( (k, self.metaFormatDict[k].getDescription() ) )
        return description_strings


    def getMetaFilterDescription(self):
        """ return a list of tuples with (id,description) of FORMAT  metalines """
        description_strings=[]
        for k in self.metaFilterDict.keys():
            description_strings.append( (k, self.metaFilterDict[k].getDescription() ) )
        return description_strings





    def setFileFormat(self, formatline):
        """ formatline is ##fileformat=VCFv4.1 """
        self.fileFormat=formatline

    def getFileFormat(self):
        return self.fileFormat

    def parseMetaInfo(self, line):
        """ parse an ##INFO line field """
        """ ##INFO<ID=ID,Number=number,Type=type,Description=description> """
        line.strip()
        infobject=InfoLine()
        pattern= '<(.*)>'
        pattern_descrip='Description=\"(.*)\"'
        
        infostring =re.search(pattern, line).groups()[0]
        
        
        if infostring == None:
            sys.stderr.write('error in parsing meta ##INFO line!\n')
            sys.exit(1)
        else:
            descriptionstring =re.search(pattern_descrip, infostring).groups()[0]
      
            infofields=infostring.split(',')
           
            for elem in infofields:
                
                if 'Description' in elem:
                    infobject.setDescription(descriptionstring)
                    break
                else:
                    (key,val)=elem.split('=')
                    if key=='ID': infobject.setId(val)
                    elif key=='Number': infobject.setNumber(val)
                    elif key=='Type': infobject.setType(val)
                    else: pass
            self.metaInfoDict[ infobject.getId() ] = infobject
        
    def parseMetaFilter(self, line):
        """ parse ##FILTER field """
        """ ##FILTER=<ID=ID,Description=description> """
        line.strip()
        filterobject=FilterLine()
        pattern='<(.*)>'
        filterstring=re.search(pattern, line).groups()[0]
        if filterstring == None:
            sys.stderr.write("error in parsing ##FORMAT line!\n")
            sys.exit(1)
        else:
            filterfields=filterstring.split(',')
            for elem in filterfields:
                (key, val)=elem.split('=',1)
                if key=='ID': filterobject.setId(val)
                elif key=='Description': filterobject.setDescription(val)
                else: pass
            self.metaFilterDict [ filterobject.getId() ] = filterobject

    def parseMetaFormat(self,line):
        """ parse ##FORMAT """
        """ ##FORMAT=<ID=ID,Number=number,Type=type,Description=description> """
        
        formatobject = FormatLine()
        pattern= '<(.*)>'
        pattern_descrip='Description=\"(.*)\"'
        formatstring =re.search(pattern, line).groups()[0]
        
        if formatstring == None:
            sys.stderr.write('error in parsing meta ##INFO line!\n')
            exit(1)
        else:
            descriptionstring =re.search(pattern_descrip, formatstring).groups()[0]
            formatfields=formatstring.split(',')
            for elem in formatfields:
                if 'Description' in elem:
                    formatobject.setDescription (descriptionstring)
                    break
                else:
                    (key,val)=elem.split('=')
                    if key=='ID': formatobject.setId(val)
                    elif key=='Number': formatobject.setNumber(val)
                    elif key=='Type': formatobject.setType(val)

                    else: pass
            self.metaFormatDict[ formatobject.getId() ] = formatobject


    def yieldPrintMetaInfoLines(self):
        """ yield string representation of meta ##INFO line objects """
        for k in self.metaInfoDict.keys():
            yield self.metaInfoDict[k].toString()
     
    def yieldPrintMetaFormatLines(self):
        """ yield string representation of meta ##FORMAT line objects """
        for k in self.metaFormatDict.keys():
            yield self.metaFormatDict[k].toString()

    def yieldPrintMetaFilterLines(self):
        """ yield string representation of meta ##FILTER line objects """
        for k in self.metaFilterDict.keys():
            yield self.metaFilterDict[k].toString()


    ######################################################################################

class HeaderLines(object):
    """       """
    def __init__(self):
        """ class to represent VCF  headerline
            #CHROM POS ID REF ALT QUAL FILTER INFO are the fixed,standard,mandatory  header columns """

        self.headercolumns=['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO'  ]
        self.hasFormatColumn=0
        self.samplelist=[]
        

    def add_format_column(self, headerline):
        """ FORMAT column is optional since not all VCFs contain genotypes; check if FORMAT is in the header line, then add it to headercolumns list """
        if 'FORMAT' in headerline:
            self.headercolumns.append('FORMAT')
            self.hasFormatColumn=1
            

    def append_samplelist(self, headerline):
        """ check to see if vcf has sample columns and if so append th list samplelist """
        
        fields=headerline.split('\t')
        
        if len(fields) > 9:
            self.samplelist=fields[9::]
        
            
    def getSampleList(self):
        return self.samplelist

    def toString(self):
        headerstring= "\t".join(self.headercolumns)
        if len ( self.samplelist ) >=1:
            samplestring="\t".join(self.samplelist)
            outstring=headerstring+"\t"+samplestring
        else:
            outstring=headerstring
        return outstring
