#!/usr/bin/env python

import libxml2, sys, os, string, optparse, math

# This file reads in expected values from the command line, and opens up a
# FrequentOutput.dat file in the outputPath directory to compare  
# if the values are within the given tolerance limits.
# To use this file from the command line, type
# python readFrequentOutput.py --help
class ReadFrequentOutput:
	def __init__( self, outputPath, valueName, expectedValue, expectedTime, toleranceTime, toleranceValue):
		self.outputPath = outputPath
		self.valueName = valueName
		self.values = []
		self.times = []
		self.timesteps = []
		self.expectedValue = string.atof(expectedValue)
		self.expectedTime = string.atof(expectedTime)
		self.toleranceTime = string.atof(toleranceTime)
		self.toleranceValue = string.atof(toleranceValue)
		self.index = -1
		self.passed = False
		self.filename = self.outputPath + "FrequentOutput.dat"

	
	def readFile( self ):
		# This function opens the FrequentOutput.dat
		# and reads the values in the file.
		if not os.path.exists( self.filename ):
			print self.filename + " file or path does not exist.\n"
			sys.exit( 1 )
			
		doc = open( self.filename )
		
		# Read the header line
		headerLine = doc.readline()
		headers = []
		# Strip out header info
		headers = headerLine.split()
		# Check that header we are interested in is included
		headerExists = False
		for i in range( len(headers) ):
			if headers[i] == self.valueName:
				headerExists = True
		if headerExists == True:		
			current = doc.readline()
	
			# Read out the values wanted out of FrequentOutput.dat to self.
			while current:
				list = []
				list = string.split(current)
				self.timesteps.append(string.atoi(list[0]))
				self.times.append(string.atof(list[1]))
				for i in range( len( headers ) - 1 ):
					if headers[i+1] == self.valueName:
						self.values.append(string.atof(list[i]))
				current = doc.readline()		
			doc.close()
		else:
			print "No header named " + self.valueName + " in " + self.filename + "\n"
			sys.exit( 1 )
	
	def compareExpectedValue( self ):
		# Search through time list to find time close to expectedTime
		for i in range(len(self.times)):
			if (math.fabs(self.times[i] - self.expectedTime) <= self.toleranceTime):
				# this is the timestep to compare.
				self.index = i
		if self.index == -1:
			self.passed = False		
		elif (math.fabs(self.values[i] - self.expectedValue) <= self.toleranceValue):
			self.passed = True
		else:
			self.passed = False		
		
		
	def writeResultFile( self ):
		# Open results file and write result line.
		tmp = ""
		if self.passed == True:
				tmp += "Value of " + self.valueName + \
				" at expected Time " + str(self.expectedTime) +" +/- " + str(self.toleranceTime) + \
				" is within tolerance " + str(self.toleranceValue) + \
				" of " + str(self.expectedValue) + "\n" 
		else:
				tmp += "Value of " + self.valueName + \
				" at expected Time " + str(self.expectedTime) +" +/- " + str(self.toleranceTime) + \
				" is not within tolerance " + str(self.toleranceValue) + \
				" of " + str(self.expectedValue) + "\n" 
				if not self.index == -1:
					tmp += " Value("+ str(self.values[self.index]) + \
					"), Time("+ str(self.times[self.index]) +\
					") at Timestep " + str(self.timesteps[self.index]) + "\n"
				else:
					tmp += "No times in range: " + str(self.expectedTime) +\
					" +/- " + str(self.toleranceTime) + "\n"
		outputFilename = self.outputPath + "/" + self.valueName + "-result.out" 						
		file = open( outputFilename, 'w' )
		file.write( tmp	)
		file.close()

	def run(self):
		# This file runs the sub-functions after self has been initialised
		self.readFile()
		self.compareExpectedValue()
		self.writeResultFile()

			
if __name__ == '__main__':
	op = optparse.OptionParser()
	op.add_option( "-o", "--outputPath", dest="outputPath", help="output Directory path" )
	op.add_option( "-n", "--valueName", dest="valueName", help="Name of value in FrequentOutput.dat file header" )
	op.add_option( "-V", "--expectedValue", dest="expectedValue", help="expected numerical Value at given time" )
	op.add_option( "-T", "--expectedTime", dest="expectedTime", help="expected Time value" )
	op.add_option( "-t", "--toleranceTime", dest="toleranceTime", help="tolerance on time value" )
	op.add_option( "-v", "--toleranceValue", dest="toleranceValue", help="tolerance on expected value at expected time" )
	(options, args) = op.parse_args()
	if not options.outputPath:
		print "outputPath required (see --help)"
		sys.exit( 1 )
	if not options.valueName:
		print "Name of value required (see --help)"
		sys.exit( 1 )
	if not options.expectedValue:
		print "expected value required (see --help)"
		sys.exit( 1 )		
	if not options.expectedTime:
		print "expected Time required (see --help)"
		sys.exit( 1 )	

	if not options.toleranceTime:
		print "tolerance on time required (see --help)"
		sys.exit( 1 )
	if not options.toleranceValue:
		print "tolerance on value required (see --help)"
		sys.exit( 1 )		
		
	checkFile=ReadFrequentOutput(options.outputPath, options.valueName, options.expectedValue, options.expectedTime, options.toleranceTime, options.toleranceValue)
	checkFile.run()
