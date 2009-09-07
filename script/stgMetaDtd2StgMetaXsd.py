#!/usr/bin/env python

import sys
import convert

def main( argv = None ):
	if argv == None:
		argv = sys.argv
	filename = argv[1]

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

if __name__=='__main__':
	main()

