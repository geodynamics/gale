#!/usr/bin/python
import createDocs
import createHTMLDocuments, createDocument
import sys, os.path, os, string, subprocess, shutil
from createDocument import createDocuments
from createDocs import createListDictionary
from createDocs import createDoxygen

def printHelpStatement(dictionary):
# specific help statement for running Doxygen ONLY
         print "Help for createDoxygen.py. This file is designed to be standalone."
         print "To run, type ./createDoxygen.py \n"
         print "To run with input options: "
         print "./createDoxygen.py DOXYGENCONFIG DOXYGENINPUTDIR DOXYGENOUTPUTDIR DOCDATAPATH\n"
         print "Defaults:"
         print "DOXYGENCONFIG: Name of Doxygen Configuration file. Default: 'Doxyfile', Current: "+ dictionary['configFile']
         print "DOXYGENINPUTDIR: Name of directory for code base. Default: '../../', Current: "+ dictionary['directoryPath']
         print "DOXYGENOUTPUTDIR: Name of directory to which Doxygen will output it's files. Default: '../../doc/', Current: "+ dictionary['docPath']
         print "DOCDATAPATH: Name of path to StGermain/doc/ directory, which has necessary file for Doxygen run. Default: "+dictionary['directoryPath']+"/StGermain/doc/', Current: "+ dictionary['docDataPath']

         print "*******WARNING********"
         print "ALL data in "+dictionary['docPath']+ " will be overwritten or removed!"
         print "**********************" 
         print "Other defaults: (Not all necessary) \n"
         print dictionary.items()


if __name__=='__main__':
    values = sys.argv
    for i in range(5):
        if len(values)< 5:
            values.append("")


    #Define input directory (relative to ~/stgUnderworldE/ )
    if values[2] != "" :
        directoryPath = os.path.realpath(values[2])
    else:
        directoryPath = os.path.realpath('../../')
    print directoryPath
    # Create the output directory. This is relative to the ./mainProj/config page.(relative to ~/stgUnderworldE/ )
    if values[3] != "":
        docPath = os.path.realpath(values[3])
    else:
        docPath = os.path.join(directoryPath, 'doc/')
    print docPath

    # createDictionary
    mainDictionary = createListDictionary("","", "", values[1], directoryPath, docPath)

    # Add other dictionary options or reset preset options
    if values[4] != "":
        mainDictionary['docDataPath'] = os.path.realpath(values[4])
    else:
        mainDictionary['docDataPath'] = os.path.realpath(os.path.join(directoryPath,'StGermain/doc/'))
    # Setup docScriptPath
    if values[4] != "":
        mainDictionary['docScriptPath'] = os.path.join(string.rstrip(os.path.realpath(values[4]), "doc/"), "script/")
    else:
        mainDictionary['docScriptPath'] = os.path.realpath(os.path.join(directoryPath,'StGermain/script/'))
    # Set up help print statement
    if ((values[1] == "help") or (values[1] == '-h') or (values[1] =='--help') or (values[1] == 'h')):
        printHelpStatement(mainDictionary)

    else:
        #check the directory now exists.
        createDocument.checkOutputPath(docPath)

        
        createDoxygen(mainDictionary)



