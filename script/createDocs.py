#!/usr/bin/env python
import createHTMLDocuments, createDocument
import sys, os.path, os, string, subprocess, shutil
from createDocument import createDocuments
from createDocument import Project
from createDocument import Meta
from createHTMLDocuments import createHTMLDocuments
from createHTMLDocuments import DivIds
## find a list of projects in the base directory, and remove obvious non-project listings
def findProjectDirectories(mainDirectory):
    directory = os.path.realpath(mainDirectory)
    # find a list of projects in the base directory
    projectNames = os.listdir(directory)
    # remove items from list that are not projects:
    # known examples: config, DOC, Doc, doc
    directoryValues = ['config', 'DOC', 'doc', 'Doc', '.hg', 'DOXYGEN', 'build', 'script', '.sconf_temp', 'scons', 'sconf_temp']
    for value in directoryValues:
        if (projectNames.count(value) > 0) :
            projectNames.remove(value)
    # Make sure only directories are included
    projList = []
    for value in projectNames:
        #print value, os.path.isdir(str(mainDirectory + '/'+value))
        if (os.path.isdir(str(mainDirectory + '/'+value)) == True):
            projList.append(value)
    print "Creating Codex pages for projects: ", projList
    return projList

## Create a dictionary that contains all data needed to create Codex and Doxygen pages
def createListDictionary(arg1, arg2, arg3, arg4, directoryPath, docPath, docDataPath):
    dictionary = {}
    # Add items to dictionary

    dictionary['directoryPath'] = os.path.realpath(directoryPath)

    dictionary['docPath'] = os.path.realpath(docPath)
    if docDataPath == "":
        dictionary['docDataPath'] = os.path.realpath(os.path.join(directoryPath,'StGermain/doc/'))
    else:
        dictionary['docDataPath'] = os.path.realpath(docDataPath)

    if docDataPath != "":
        dictionary['docScriptPath'] = os.path.join(string.rstrip(os.path.realpath(dictionary['docDataPath']), "doc/"), "script/")
    else:
        dictionary['docScriptPath'] = os.path.realpath(os.path.join(directoryPath,'StGermain/script/'))

    # Define Codex values
    # Codex Subdir
    dictionary['codexSubDir'] = 'Codex'
    dictionary['codexSubDirPath'] = os.path.join(dictionary['docPath'], dictionary['codexSubDir'])    
    # Define stylesheet data to input
    dictionary['stylesheetList'] = [os.path.join(dictionary['docDataPath'], 'print.css'),os.path.join(dictionary['docDataPath'], 'menu.css'),os.path.join(dictionary['docDataPath'], 'codexStylesheet.css')]

    # Define javascript to input
    dictionary['script'] =[os.path.join(dictionary['docDataPath'], 'menuscript.js')] 
    #Define blurb
    dictionary['documentBlurb'] ="This is a list of the available components."
    
    # Define picture directory
    dictionary['imagePath'] = os.path.join(dictionary['docDataPath'], 'images')

    #Define output picture directory
    dictionary['pictureDirectory'] = 'images/'

    dictionary['arg1'] = arg1
    dictionary['arg2'] = arg2
    dictionary['arg3'] = arg3
    dictionary['arg4'] = arg4

    # Get data read in from command line
        
    if (dictionary['arg1'] == ""):
        dictionary['metaFlag'] = 'dtd'
    elif ((dictionary['arg1'] == "help") or (dictionary['arg1'] == '-h') or (dictionary['arg1'] =='--help') or (dictionary['arg1'] == 'h')):
        dictionary['metaFlag'] = 'dtd'
    else:
        dictionary['metaFlag'] = arg1

    if (dictionary['arg2'] == ""):
        dictionary['extWeb'] = [['http://www.underworldproject.org','Underworld Home Page']]
    else:
        dictionary['extWeb'] = arg2
    
    if (dictionary['arg3'] == ""):
        dictionary['email'] = 'underworld-users@vpac.org'
    else:
        dictionary['email'] = arg3


    if (dictionary['arg4'] == ""):
            dictionary['configFile'] = 'Doxyfile'
    else:
            dictionary['configFile'] = str(arg4)
    
        #Define location of content for index html page
    if (dictionary['metaFlag'] == 'dtd'):
        dictionary['indexFileName'] = os.path.join(dictionary['docDataPath'], 'IndexContent-Dtd.html')
    elif (dictionary['metaFlag'] == 'xsd'):
        dictionary['indexFileName'] = os.path.join(dictionary['docDataPath'], 'IndexContent-Xsd.html')
    else:
        dictionary['indexFileName'] = ""

    # Doxygen values
    dictionary['doxygenSubDir'] = 'Doxygen'
    dictionary['doxygenSubDirPath'] = os.path.join(dictionary['docPath'], dictionary['doxygenSubDir'])
    dictionary['configPath'] = os.path.join(dictionary['docDataPath'], dictionary['configFile'])
    dictionary['headerFilterPath'] = dictionary['docScriptPath']+"/headerfilter.py"
    # All the following values are contained within 'docDataPath':
    dictionary['htmlHeader'] = "header.html"
    dictionary['htmlFooter'] =  "footer.html"
    dictionary['htmlStylesheet'] =  "customdoxygen.css"
    dictionary['htmlImagesPath'] = "doxyimage/"
    dictionary['projectNumber'] = 'Bleeding Edge'
    dictionary['configNew'] =  "Doxyfile.new"
    dictionary['htmlDir'] = "html/"
    return dictionary

## Print help statement appropriate for 'createDocs.py'
def printHelpStatement(dictionary):
         print "Help for createDocs.py. This file is designed to be run with scons doc."
         print "Use createHTMLDocs.py for the standalone version. "
         print "To run, type ./createDocs.py \n"
         print "To run with input options: "
         print "./createDocs.py METATYPE WEBPAGES EMAIL DOXYGENCONFIG\n"
         print "Defaults:"
         print "METATYPE: 'dtd' or 'xsd'"
         print "WEBPAGES: [['http://www.underworldproject.org','Underworld Home Page']]"
         print "EMAIL: 'underworld-users@vpac.org'\n"
         print "DOXYGENCONFIG: Doxyfile"

         print "*******WARNING********"
         print "ALL data in "+dictionary['docPath']+ " will be overwritten or removed!"
         print "**********************" 
         print "Other defaults: \n"
         print dictionary.items()

## Main file to run component codex creator with defaults built in
def createCodex(dictionary):

    #Now create the projects.
    
    projectList = []

    projectNames = findProjectDirectories(dictionary['directoryPath'])

    for projectName in projectNames:
        # check project name in current path already
        if (string.count(dictionary['directoryPath'],  str("/"+projectName)) > 0) :
           path =  os.path.realpath(dictionary['directoryPath']) + "/"
        # check if adding projectName
        else:
           path = os.path.realpath(dictionary['directoryPath']) + "/"+ str(projectName) + "/"
        print "Creating data for " + projectName
        project =Project(projectName, path,  dictionary['metaFlag'])
        project.assignMetas()
        projectList.append(project)
        projectList.sort()
        


    #now create the HTML documents
    print "Now creating HTML documents in directory: " + os.path.realpath(dictionary['codexSubDirPath'])
    
    htmlDocuments = createHTMLDocuments("Component Codex", dictionary['documentBlurb'],   projectList, dictionary['codexSubDirPath'],  dictionary['stylesheetList'],  dictionary['extWeb'],  dictionary['email'], dictionary['script'],  dictionary['pictureDirectory'],  dictionary['indexFileName'],  "Google", "")

    htmlDocuments.createHTMLPages()
    # These two aren't necc. for this file, as the images are already in the right place.    
    htmlDocuments.copyPictures(dictionary['imagePath'])
    htmlDocuments.copyStylesheets()

## create the doxygen web pages
def createDoxygen(dictionary):

    #Todo: check doxygen exists

 
    #Todo: Create directory to write Doxygen to
    if os.path.isdir(dictionary['doxygenSubDirPath']) == False:
        print "Path for Doxygen output does not exist. Creating directory: " + dictionary['doxygenSubDirPath'] 
        os.mkdir(dictionary['doxygenSubDirPath'])
    #Check whether Doxyfile exists in dictionary[docPath]

    if os.path.isfile(dictionary['configPath'])==True:
        print "Config file, "+dictionary['configPath']+ " exists."
    elif os.path.isfile(dictionary['configPath'])==False:
        print "Config file " + str(dictionary['configPath']) + " does not exist."
        print "Creating Doxygen config file, name: " + dictionary['configFile']

    # Modify doxygen entries to work with current directory system
   
    # Open Doxyfile
    doxyfile = open(dictionary['configPath'])   
    doxyfileString = doxyfile.readlines()
    doxyfile.close()
    newDoxyfileString = ""
    for line in doxyfileString:
        newLine = ""
        
        # Find the "INPUT" field and add the root project dir to it
        if (line.find("INPUT=") != -1) or ((line.find("INPUT ") != -1) and (line.find("=") != -1)):
            newLine = "INPUT = "+ dictionary['directoryPath']
            print "Found INPUT"
            print newLine
        elif (line.find("STRIP_FROM_PATH") != -1)and (line.find("=") != -1):
            newLine = "STRIP_FROM_PATH = "+ dictionary['directoryPath']
            print "Found STRIP_FROM_PATH"
            print newLine
        elif ((line.find("PROJECT_NUMBER") != -1) and (line.find("=") != -1)):
            newLine = "PROJECT_NUMBER = "+ dictionary['projectNumber']
            print "Found PROJECT_NUMBER"
            print newLine
        elif ((line.find("OUTPUT_DIRECTORY") != -1) and (line.find("=") != -1)):
            newLine = "OUTPUT_DIRECTORY = " + dictionary['docPath'] +"/"+dictionary['doxygenSubDir']+'/'
            print "Found OUTPUT_DIRECTORY"
            print newLine
        elif ((line.find("FILTER_PATTERNS") != -1) and (line.find("=") != -1)):
            newLine = "FILTER_PATTERNS = *.h="+dictionary['headerFilterPath']
            print "Found FILTER_PATTERNS: "
            print newLine
        elif ((line.find("HTML_STYLESHEET") != -1) and (line.find("=") != -1)):
            newLine = "HTML_STYLESHEET = "+dictionary['docDataPath']+"/"+dictionary['htmlStylesheet']
            print "Found HTML_STYLESHEET: "
            print newLine
        elif ((line.find("HTML_FOOTER") != -1) and (line.find("=") != -1)):
            newLine = "HTML_FOOTER = "+dictionary['docDataPath']+"/"+dictionary['htmlFooter']
            print "Found HTML_FOOTER: "
            print newLine
        elif ((line.find("HTML_HEADER") != -1) and (line.find("=") != -1)):
            newLine = "HTML_HEADER = "+dictionary['docDataPath']+"/"+dictionary['htmlHeader']
            print "Found HTML_HEADER: "
            print newLine

        else:
            newLine = line
 
        newDoxyfileString = newDoxyfileString + newLine

    
    #write new doxyfile
    
    doxyfileNew = open(os.path.join(dictionary['docDataPath'],dictionary['configNew']), 'w')
    doxyfileNew.write(newDoxyfileString)
    doxyfileNew.close()
    
    # copy pics to new directory
    copyPics(os.path.join(dictionary['docDataPath'],dictionary['htmlImagesPath']), os.path.join(dictionary['doxygenSubDirPath'],dictionary['htmlDir']))    

    #run doxygen ( this command checks PATH for doxygen.
    print "Creating Doxygen documentation for " + str(dictionary['configPath'])
    print " using modified file: "+str(os.path.join(dictionary['docDataPath'],dictionary['configNew'])) 
    print "Output located in " + dictionary['doxygenSubDirPath'] 
    os.execlp('doxygen',  'doxygen', str(os.path.join(dictionary['docDataPath'],dictionary['configNew'])) )



    #remove temporary Doxyfile
    os.remove(str(os.path.join(dictionary['docDataPath'],dictionary['configNew'])))

## Basic copy function for directory contents. Is not recursive.
def copyPics(pathFrom, pathTo):
    print "Copying content of directory of pictures :" + os.path.realpath(pathFrom) + "/* to directory: " + os.path.realpath(pathTo) 
    createDocument.checkOutputPath(pathTo)
    # Copy everything there to the moment
    os.system('cp '+ os.path.realpath(pathFrom) + "/* "+os.path.realpath(pathTo) )
    
   
if __name__=='__main__':
    values = sys.argv
    for i in range(5):
        if len(values)< 5:
            values.append("")
    

    #Define input directory (relative to ~/stgUnderworldE/ )
    directoryPath = os.path.realpath('.')
    print directoryPath
    # Create the output directory. This is relative to the ./mainProj/config page.(relative to ~/stgUnderworldE/ )
    docPath = os.path.realpath('./doc/')
    
    # createDictionary
    mainDictionary = createListDictionary(values[1],values[2], values[3], values[4], directoryPath, docPath, "")

    # Set up help print statement
    if ((mainDictionary['arg1'] == "help") or (mainDictionary['arg1'] == '-h') or (mainDictionary['arg1'] =='--help') or (mainDictionary['arg1'] == 'h')):
        printHelpStatement(mainDictionary)

    else:
        #check the directory now exists.
        createDocument.checkOutputPath(docPath)

        createCodex(mainDictionary)
        createDoxygen(mainDictionary)


