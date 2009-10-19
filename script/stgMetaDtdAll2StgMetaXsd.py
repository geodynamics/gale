#!/usr/bin/env python

import os, sys
import convert
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

			xsd = convert.dtd2MetaXsd( xml_text )
			try:
				print( xsd )
			except UnicodeEncodeError:
				print 'File contains unicode characters... printing as c-string'
				print( xsd.encode( 'utf-8' ) )

		except:
			failed.append( filename )
	
	if len( failed ):
		print 'Failed list (' + str( len( failed ) ) + ' of ' + str( len( filenameList ) ) + ')... ' 
		pprint.pprint( failed )
	else:
		print 'All ' + str( len( filenameList ) ) + ' passed conversion test'

if __name__=='__main__':
	main()

