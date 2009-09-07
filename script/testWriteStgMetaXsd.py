#!/usr/bin/env python

import stgMetaXsd
import testData
import pprint

def main():
	d = testData.stgMetaXsdDict()
	doc = stgMetaXsd.createXML( d )

	s = doc.toprettyxml()
	print s
	#pprint.pprint( s ) # used to get text for testData.stgMetaXsdXml

	if s == testData.stgMetaXsdXml():
		print 'Passed'
	else:
		print 'Failed'


if __name__=='__main__':
	main()

