#!/usr/bin/python
import os, os.path, sys,  datetime
import convert,  util,  stgDtd,  stgMetaXsd

## checkOutputPath checks that the value given is a valid path.    
def checkOutputPath(directory):
    #check if path is valid
    if (os.path.isdir(directory)==False):
        print(directory+" does not exist ... Making directory")
        os.mkdir(directory)


## createDocuments is a generic creation class designed to be inherited
# by the specific document creation classes. ie createHTMLDocuments or
# createPDFDocuments etc.
class createDocuments():
    def __init__(self,  documentName,  projectList, documentPath):
        self.name= documentName
        self.projectList = projectList
        self.documentList =[]
        self.path = os.path.realpath(documentPath)
        self.indexPage = []

## The Project class is an all-inclusive class for all the project content.
# This includes the list of associated meta files, it's path, name, metatype etc.
class Project():
    def __init__(self,  name, path,  metaFlag):
        self.name = name
        self.metas = []
        self.path =path
        #string = either 'xsd' or 'dtd'
        self.metaFlag = metaFlag
        
        #check if path is valid
        if (os.path.isdir(self.path)==False):
            sys.exit(self.path + " is not a valid directory.")
    ## This function assigns the metas to the project. It opens all meta files
    # and creates a dictionary based on the metaType specified.
    def assignMetas(self):
        #search path for all .meta files
        #TODO check if this accesses the .hg or .svn dirs. If so, might have to prune. may use code in checkspelling.py
        fileList = util.dirEntries( self.path, True, 'meta' )
        
        #open each meta file, extract out text and convert to meta dictionary
        for filepath in  fileList:
            if os.path.isfile(filepath):
                file =open(filepath)
                text = file.read()

                metaFile = Meta(text)
                if self.metaFlag == 'xsd':
                    metaFile.createMetaDictionaryXsd()
                    # add meta dictionary to list if it is part of project
                    if metaFile.dictionary['info']['subject'] == self.name:
                        self.metas.append(metaFile)
                elif self.metaFlag =='dtd': 
                    metaFile.createMetaDictionaryDtd()
                    # add meta dictionary to list if it is part of project
                    if metaFile.dictionary['Project'] == self.name:
                        self.metas.append(metaFile)                   
        self.metas.sort()
## The Meta Class is just a glorified dictionary with some functions.
class Meta():
    def __init__(self,  xmlText):
        self.dictionary = {}
        self.xmlText = str(xmlText)
    ## Create the classic dictionary type.    
    def createMetaDictionaryDtd(self):
        
        # Parse DTD
        try:
            self.dictionary = stgDtd.readXML( self.xmlText )
        except:
            print 'Failed to parse as a StGermain DTD'
            raise
        
    ## create the new dictionary type.
    def createMetaDictionaryXsd(self):
        self.createMetaDictionaryDtd()
        # Convert DTD-dict to XSD-dict
        try:
            self.dictionary = convert.dtdDict2metaXsdDict( self.dictionary )
        except:
            print 'Failed to convert information from a StGermain Meta DTD to a StGermain Meta XSD'
            raise
            
            

