#!/usr/bin/env python

import convert
import testData
import pprint

def main():
	dtd = testData.stgMetaDtdDict()
	xsd = convert.dtdDict2metaXsdDict( dtd )

	pprint.pprint( xsd )

	test = testData.stgMetaXsdDict()
	pprint.pprint( test )

	if test == xsd:
		print 'Passed'
	else:
		print 'Failed'


if __name__=='__main__':
	main()

