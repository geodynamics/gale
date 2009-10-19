#!/usr/bin/env python

import os, sys
import pprint
import stgDtd
import util

def main( argv = None ):
	if argv == None:
		argv = sys.argv
	dirname = argv[1]

	failed = []
	filenameList = util.dirEntries( dirname, True, 'meta' )
	for filename in filenameList:
		print( filename + '...' )
		try:
			# Read the file into a string
			xml_file = file( filename )
			xml_lines = xml_file.readlines()
			xml_text = ""
			for l in xml_lines:
				xml_text += str(l)

			dtdDict = stgDtd.readXML( xml_text )
			pprint.pprint( dtdDict )
		except:
			failed.append( filename )
	
	if len( failed ):
		print 'Failed list (' + str( len( failed ) ) + ' of ' + str( len( filenameList ) ) + ')... ' 
		pprint.pprint( failed )
	else:
		print 'All ' + str( len( filenameList ) ) + ' passed parsing test'

if __name__=='__main__':
	main()

