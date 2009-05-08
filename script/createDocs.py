#!/usr/bin/python
import createHTMLDocuments, createDocument
import sys, os.path, os, string, subprocess
from createDocument import createDocuments
from createDocument import Project
from createDocument import Meta
from createHTMLDocuments import createHTMLDocuments
from createHTMLDocuments import DivIds
def findProjectDirectories(mainDirectory):
    directory = os.path.realpath(mainDirectory)
    # find a list of projects in the base directory
    projectNames = os.listdir(directory)
    # remove items from list that are not projects:
    # known examples: config, DOC, Doc, doc
    directoryValues = ['config', 'DOC', 'doc', 'Doc', '.hg', 'DOXYGEN', 'build']
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


def createListDictionary(arg1, arg2, arg3, arg4, directoryPath, docPath):
    dictionary = {}
    # Add items to dictionary

    dictionary['directoryPath'] = directoryPath

    dictionary['docPath'] = docPath

    # Define Codex values
    # Codex Subdir
    dictionary['codexSubDir'] = 'Codex'
    dictionary['codexSubDirPath'] = os.path.join(dictionary['docPath'], dictionary['codexSubDir'])    
    # Define stylesheet data to input
    dictionary['stylesheetList'] = [os.path.join(docPath, 'print.css'),os.path.join(dictionary['docPath'], 'menu.css'),os.path.join(dictionary['docPath'], 'codexStylesheet.css')]

    # Define javascript to input
    dictionary['script'] =[os.path.join(dictionary['docPath'], 'menuscript.js')] 
    #Define blurb
    dictionary['documentBlurb'] ="This is a list of the available components."
    
    # Define picture directory
    dictionary['imagePath'] = os.path.join(dictionary['docPath'], 'images')

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
        dictionary['indexFileName'] = os.path.join(dictionary['docPath'], 'IndexContent-Dtd.html')
    elif (dictionary['metaFlag'] == 'xsd'):
        dictionary['indexFileName'] = os.path.join(dictionary['docPath'], 'IndexContent-Xsd.html')
    else:
        dictionary['indexFileName'] = ""

    # Doxygen values
    dictionary['doxygenSubDir'] = 'Doxygen'
    dictionary['doxygenSubDirPath'] = os.path.join(dictionary['docPath'], dictionary['doxygenSubDir'])
    dictionary['configPath'] = os.path.join(dictionary['docPath'], dictionary['configFile'])
    dictionary['headerFilterPath'] = dictionary['docPath']+"/headerfilter.py"
    dictionary['projectNumber'] = 'Bleeding Edge'
    dictionary['configPathNew'] = os.path.join(dictionary['docPath'], "Doxyfile.new")
    return dictionary

# Print help statement
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

         print "Other defaults: \n"
         print dictionary.items()

# Main file to run component codex creator with defaults built in
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
        elif (line.find("STRIP_FROM_PATH") != -1)and (line.find("=") != -1):
            newLine = "STRIP_FROM_PATH = "+ dictionary['directoryPath']
            print "Found STRIP_FROM_PATH"
        elif ((line.find("PROJECT_NUMBER") != -1) and (line.find("=") != -1)):
            newLine = "PROJECT_NUMBER = "+ dictionary['projectNumber']
            print "Found PROJECT_NUMBER"
        elif ((line.find("OUTPUT_DIRECTORY") != -1) and (line.find("=") != -1)):
            newLine = "OUTPUT_DIRECTORY = " + dictionary['docPath'] +"/"+dictionary['doxygenSubDir']+'/'
            print "Found OUTPUT_DIRECTORY"
            print newLine
        elif ((line.find("FILTER_PATTERNS") != -1) and (line.find("=") != -1)):
            newline = "FILTER_PATTERNS = *.h="+dictionary['headerFilterPath']
            print "Found FILTER_PATTERNS"

        #HTML_STYLESHEET        =
        #HTML_FOOTER            =
        #HTML_HEADER            =

        else:
            newLine = line
 
        newDoxyfileString = newDoxyfileString + newLine

    
    #write new doxyfile

    doxyfileNew = open(dictionary['configPathNew'], 'w')
    doxyfileNew.write(newDoxyfileString)
    doxyfileNew.close()

    #run doxygen ( this command checks PATH for doxygen.
    print "Creating Doxygen documentation for " + str(dictionary['configPath'])
    print " using modified file: "+str(dictionary['configPathNew']) 
    print "Output located in " + dictionary['doxygenSubDirPath'] 
    os.execlp('doxygen',  'doxygen', str(dictionary['configPathNew']) )

    #remove temporary Doxyfile
    os.remove(dictionary['configPathNew'])
    
#def createIndexPage(dictionary):

   #
if __name__=='__main__':
    values = sys.argv
    for i in range(5):
        if len(values)< 5:
            values.append("")
    
    #print "values =" + str(values)

    #Define input directory (relative to ~/stgUnderworldE/ )
    directoryPath = os.path.realpath('.')
    print directoryPath
    # Create the output directory. This is relative to the ./mainProj/config page.(relative to ~/stgUnderworldE/ )
    docPath = os.path.realpath('./doc/')
    createDocument.checkOutputPath(docPath)

    # createDictionary
    mainDictionary = createListDictionary(values[1],values[2], values[3], values[4], directoryPath, docPath)

    # Set up help print statement
    if ((mainDictionary['arg1'] == "help") or (mainDictionary['arg1'] == '-h') or (mainDictionary['arg1'] =='--help') or (mainDictionary['arg1'] == 'h')):
        printHelpStatement(mainDictionary)

    else:
        createCodex(mainDictionary)
        createDoxygen(mainDictionary)
#        createIndexPage(mainDictionary)

