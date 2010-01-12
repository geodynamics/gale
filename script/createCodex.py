#!/usr/bin/python

## This is a codelet to run just the codex with standard defaults from anywhere in the code.
import createDocs
import createHTMLDocuments, createDocument
import sys, os.path, os, string, subprocess, shutil
from createDocument import createDocuments
from createDocs import createListDictionary
from createDocs import createCodex

## Function for specific help statement for running Doxygen ONLY
def printHelpStatement(dictionary):

         print "Help for createDoxygen.py. This file is designed to be standalone."
         print "To run, type ./createCodex.py \n"
         print "To run with input options: "
         print "./createCodex.py  METATYPE WEBPAGES EMAIL CODEXINPUTDIR CODEXOUTPUTDIR \n"
         print "Defaults:"
         print "METATYPE: 'dtd' or 'xsd'"
         print "WEBPAGES: [['http://www.underworldproject.org','Underworld Home Page']]"
         print "EMAIL: 'underworld-users@vpac.org'\n"
         print "CODEXINPUTDIR: Name of base directory.  Default: '../../', Current: "+ dictionary['directoryPath']
         print "CODEXOUTPUTDIR: Name of directory to which Codex will output it's files. Default: '../../doc/', Current: "+ dictionary['docPath']
         print "DOCDATAPATH: Name of path to StGermain/doc/ directory, which has necessary files for the Codex to run. Default: "+dictionary['directoryPath']+"/StGermain/doc/', Current: "+ dictionary['docDataPath']
         print "*******WARNING********"
         print "This codex file ONLY works with Stg base directory. "
         print " To only create codex for smaller samples, use createHTMLDocuments.py"
         print "**********************" 
         print "*******WARNING********"
         print "ALL data in "+dictionary['docPath']+ " will be overwritten or removed!"
         print "**********************" 
         print "Other defaults: (Not all necessary) \n"
         print dictionary.items()

if __name__=='__main__':
    values = sys.argv
    for i in range(7):
        if len(values)< 7:
            values.append("")
    

    
    #Define input directory (relative to ~/stgUnderworldE/ )
    if values[4] != "" :
        directoryPath = os.path.realpath(values[4])
    else:
        directoryPath = os.path.realpath('../../')
    print directoryPath

    # Create the output directory. This is relative to the ./mainProj/config page.(relative to ~/stgUnderworldE/ )
    if values[5] != "":
        docPath = os.path.realpath(values[5])
    else:
        docPath = os.path.join(directoryPath, 'doc/')
    print docPath

    # createDictionary
    mainDictionary = createListDictionary(values[1],values[2], values[3], "", directoryPath, docPath, values[6])

    # Set up help print statement
    if ((mainDictionary['arg1'] == "help") or (mainDictionary['arg1'] == '-h') or (mainDictionary['arg1'] =='--help') or (mainDictionary['arg1'] == 'h')):
        printHelpStatement(mainDictionary)

    else:
        #check the directory now exists.
        createDocument.checkOutputPath(docPath)

        createCodex(mainDictionary)

