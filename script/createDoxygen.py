#!/usr/bin/python
# The following documentation is designed to sit at the top of the main doxygen page!
## \mainpage 
## Documentation in these Doxygen pages is designed for developers.
## There are the following sections: Classes, Files and Directories
## \section Classes
## Extracted from C and python documentation, this has details on all classes defined in
## the code, including 'pseudo-classes' based on macros.
## \section Files
## Contains the extracted *.c, *.h and *.py files.
## \section Directories
## Contains a listing of all directories in the code base.
## \section Namespaces
## This is where functions
## not contained in a class but referenced by file name will appear.
## \section Related Pages
## Contains the Todo page - a list of extracted Todo's from all the code base; 
## and The Deprecated List, which contains deprecated classes, functions and files. 

import createDocs
import createHTMLDocuments, createDocument
import sys, os.path, os, string, subprocess, shutil
from createDocument import createDocuments
from createDocs import createListDictionary
from createDocs import createDoxygen

## Function for specific help statement for running Doxygen ONLY
def printHelpStatement(dictionary):

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
    mainDictionary = createListDictionary("","", "", values[1], directoryPath, docPath, values[4])

    # Set up help print statement
    if ((values[1] == "help") or (values[1] == '-h') or (values[1] =='--help') or (values[1] == 'h')):
        printHelpStatement(mainDictionary)

    else:
        #check the directory now exists.
        createDocument.checkOutputPath(docPath)

        
        createDoxygen(mainDictionary)



