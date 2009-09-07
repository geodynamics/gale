#!/usr/bin/env python

import xml.sax
import stgDtd
import pprint
import sys

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

	dtdDict = stgDtd.readXML( xml_text )
        pprint.pprint( dtdDict )

if __name__=='__main__':
	main()

