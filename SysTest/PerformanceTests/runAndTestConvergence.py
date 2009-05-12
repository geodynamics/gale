#! /usr/bin/env python

import sys
import re
import os
import commands

######  Main program    ######

# 1) Run the same xml at the different resolutions
#    these will generate error data

# 2) Use a linear regression algorithm on the data generated
#    and check the results for good behavior and acceptable convergence

######  End main program ######

def runTests():

	# read commandLine args
	ii = 0	
	xmlFile = "NULL"	
	optFile = "NULL"
	testsToRun = 4
	procs = [1,1,1,1]
	commandLines = ["--elementResI=16 --elementResJ=16 --elementResK=16","--elementResI=32 --elementResJ=32 --elementResK=32",
"--elementResI=64 --elementResJ=64 --elementResK=64",
"--elementResI=128 --elementResJ=128 --elementResK=128"]

	# check if xml exists and options file is specified 	
	for arg in sys.argv:
		RE1 = re.compile( '.*\.xml$' )
		RE2 = re.compile( '^-optionsFile*' )

		RESULT1 = RE1.match( arg )
		RESULT2 = RE2.match( arg )

		if RESULT1 != None:
			xmlFile = arg
		elif RESULT2 != None:
			optFile = arg

	if xmlFile == "NULL":
		sys.stdout.write( "\nNo xml file specified, stopped" )
		sys.exit(0)
	
	if os.path.exists( xmlFile ) != True:
		sys.stdout.write( "\nCannot find input file: " + xmlFile + ", stopped" )
		sys.exit(0)

	# check if options file is given, otherwise run default
	if optFile == "NULL":
		sys.stdout.write( "\nNo run options file specifed using commandline arg '-optionsFile xxx'\nUsing default serial mode, with commandLines:\n" )
	
		for comm in commandLines:
			sys.stdout.write( comm + "\n" ) 
	else:
		if os.path.exists( optFile ) != True:
			sys.stdout.write( "\nCannot find run options file $optFile, stopped" )
			sys.exit(0)

	# do checks on executables and log files
	executable = "udw"
	# Need to check for an executable
	if os.path.exists( "../../../build/bin/StGermain" ) != True:
		sys.stdout.write( "\nCan't find ../../../build/bin/StGermain - the executable which runs the test, stopped" )
		sys.exit(0)

	sys.stdout.write( "\n--- Testing the convergence of " + xmlFile + "---\n" )
	# is the symbolic link there, if not create it
	if os.path.exists( executable ) != True:
		command = "ln -s ../../../build/bin/StGermain " + executable
		sys.stdout.write( "\n" + command + "\n\n" )
		os.system( command )

	# check if there's a log dir
	if os.path.exists( "log/" ) != True:
		command = "mkdir log"
		os.system( command )

	stdout = "log/" + xmlFile + "_runs.stdout"
	stderr = "log/" + xmlFile + "_runs.stderr"

	# remove old log file, if it exists
	if os.path.exists( stdout ) == True:
		command = "rm " + stdout
		os.system( command ) 

	# remove old cvg file, if it exists
	command = "ls *\.cvg 2>/dev/null"
	cvg = os.system( command )

	if cvg == 0:
		command = "rm *.cvg"
		os.system( command )
	
	#commence running
	run_I = 0

	for run_I in range (0, testsToRun):
		# run test case
		command = "mpirun -np " + str( procs[run_I] ) + " ./" + executable + " " + xmlFile + " " + commandLines[run_I] + " --pluginData.appendToAnalysisFile=True >" + stdout
		command = command + " 2>" + stderr
		os.system( command )
		sys.stdout.write( command )
		# check error stream for error result

		try:
			ERROR = open( stderr, "r" )
		
			try:	
				line = ERROR.readline()

				while line:
					RE = re.compile( '[E][e]rror' )
					RESULT = RE.match( line )

					if RESULT != None:
						ERROR.close()
						sys.stdout.write( "\n\tError: see " + stderr + " or " + stdout + " - stopped" )
						os.exit(0)
 	
					line = ERROR.readline()	
			finally:
				ERROR.close()
				command = "rm " + stderr + "; " + command
				sys.stdout.write( "....finished\n" )

		except IOError:
			pass

	# removing softlink
	command = "rm " + executable
	sys.stdout.write( command )
	os.system( command )
	sys.stdout.write( "\n\n--- Finished convergence runs ---\n\n" )


def readOptionsFile( optFile_param, procs_param, commandLines_param ):
	# $line_I represents the number of tests to run	
	line_I = 0	
	# open options file

	try:
		OPTFILE = open( optFile_param, "r" )
	
		try:
			line = OPTFILE.readline()

			while line:
				line = OPTFILE.readline()
				# only process lines that start with np
		finally:
			OPTFILE.close()
	except IOError:
		pass
		
def testNumbersAgainstExpected( datFile, expectedFile ):
	ii = 0
	needLabels = 1
	jj = 0
	input = list()
	keys = list()
	expected = list()

	try:	
		INPUT = open( datFile, "r" )
	
		try:
			sys.stdout.write( "\n--- Testing the output against the expected file " + expectedFile + " ---\n\n" )

			# get results from inputfile
			line = INPUT.readline()

			while line:
				RE = re.compile( '^\-\-\-\-\-\-\-\-\-' )
				RESULT = RE.match( line )

				if RESULT != None:
					break
	
				RE = re.compile( '^#Res\s+(.*)' ) 
				RESULT = RE.match( line )

				if (RESULT != None) and (needLabels == 1):
					keys = line.split( '/\s+/' )	
					needLabels = 0	

				RE = re.compile( '^\s*\D' )
				RESULT = RE.match( line )

				if RESULT != None:
					continue

				RE = re.compile( '^\s*\d' )
				RESULT = RE.match( line )

				if RESULT != None:
					input[ii] = [ line.split( '/\s/' ) ];
					ii = ii + 1
		finally:
			INPUT.close()
	except IOerror:
		pass

	# get results from expected file
	ii = 0

	try:
		EXPECTED = open( expectedFile, "r" )
	
		try:
			line = EXPECTED.readline()

			while line:
				RE = re.compile( '^\-\-\-\-\-\-\-\-\-' )
				RESULT = RE.match( line )

				# exit input if ------ is detected
				if RESULT != None:
					break
		
				RE = re.compile( '^\s*\D' )
				RESULT = RE.match( line )

				# disregard lines starting with non-Digits	
				if RESULT != None:
					continue

				RE = re.comile( '^\s*\d' )
				RESULT = RE.match( line )

				if RESULT != None:
					expected[ii] = [ line.split( '/\sa/' ) ]
					ii = ii + 1

		finally:
			EXPECTED.close()
	except IOError:
		pass
		
	for index, item in enumerate( input ):
		if item != expected[ index ]:
			sys.stdout.write( "Error: The expected file \n\"" + expectedFile + "\nhas a different number of results than the test file \n\"" + datFile + "\nso the results testing can't be done.\nRemove argument '-againstExpected' to not check the resulting errors agains the expected file.\n" )
			os.exit(0)

	error = 0.0
	r_error = 0.0
	r_tol = 0.05

	sys.stdout.write( "relative tolerance is set to " + r_tol + "%\n" )

	for ii in range( 0, len(input) ):
		for jj in range( 1, len(keys) ):
			error = abs( input[ii][jj] - expected[ii][jj] )
			r_error = 100 * ( error / expected[ii][jj] )

			if r_error > r_tol:
				sys.stdout.write( "Results of " + keys[jj-1] + " differs by " + r_error + "\n" )

			# This is to verbose be could be useful to understand error
			else:
				sys.stdout.write( "Results of " + keys[jj-1] + " are " + r_error )
			
	sys.stdout.write( "All values within tolerance\n" )
	sys.stdout.write( "\n-- Finished testing the output against the expected file ---\n" )
	
def generateConvergence():
	# check if the cvg can be found
	command = "ls *\.cvg 2>/dev/null"
	datFile = ""
	cvgStatus = os.system( command )

	if cvgStatus != 0:
		sys.stdout.write( "There is no cvg file found analyse\n" )
		sys.exit(0)

	cvgFile = commands.getoutput( "ls *.cvg" )	
	# run David A. May's convergence tester, cheers Mayhem		
	command = "./generate-convergence-rates.pl -errorfile " + cvgFile + " -graphs"
	os.system( command )

	datFile = cvgFile + "-convergence-rates.dat"
	if os.path.exists( datFile ) != True:
		sys.stdout.write( "Can't find the convergence rate file:" + datFile )
		sys.exit(0)

	runAgainstExpected = 0

	for arg in sys.argv:
		RE = re.compile ( '-againstExpected' )
		RESULT = RE.match( arg )

		if RESULT != None:
			runAgainstExpected = 1
			break

	if runAgainstExpected == 1:
		if os.path.exists( "./expected/" + datFile ) == True:
			testnumbersAgainstExpected( datFile, "./expected/" + datFile )
		else:
			sys.stdout.write( "Can't find the expected file: ./expected/" + datFile + "\n Consider adding argument \'-convergeOnly\' to the executable to not check the resulting errors agains the expected file.\n" )
			sys.exit(0)

	# remove the old file
	command = "rm *.cvg"
	#os.system( command )

	return datFile

def testConvergence( datFile_param ):
	command = ""
	keys = list()
	cvgRates = list()
	correlations = list()
	container = list()
	needLabels = 1
	line = ""

	#test convergence numbers

	try:
		INPUT = open( datFile_param, "r" )
	
		try:
			line = INPUT.readline()

			while line:
				line = line.rstrip()
		
				# parse for convergence rates
				RE = re.compile( '^cvg\.\srate\s+(.*)' )
				RESULT = RE.match( line )

				if RESULT != None:
					for object in line.split():
						if object != " ":
							cvgRates.append( object )
					#cvgRates = line.split( '\s+' )

				# parse for correlation results
				RE = re.compile( 'corr' )
				RESULT = RE.match( line )

				if RESULT != None:
					for object in line.split():
						if object != " ":
							correlations.append( object )
				#correlations = line.split( " " )

				# parse for variable labels
				RE = re.compile( '^#Res\s+(.*)' )
				RESULT = RE.match( line )

				if ( RESULT != None ) and ( needLabels == 1 ):
					for object in line.split():
						if object != " ":
							keys.append( object )	
							needLabels = 0	
				#keys = line.split( '\s+' )

				line = INPUT.readline()

			# now check results and report
			ii = 0
			jj = 0
			length = 0
			status = 0
			result = "Pass"
			report = ""
			padding = ""

			sys.stdout.write( "\n--- Results of convergence ---\n\n" )
	
			for ii in range( 1, len( correlations ) ):
				status = "good"
				report = ""
				padding = ""

				length = 25-len( keys[ii] )

				for jj in range( 0, length ):
					padding += " "

				sys.stdout.write( keys[ii] + " " + padding + " - convergence rate = " + cvgRates[ii+1] + " - correlation" + correlations[ii] + " .... " )

				if correlations[ii] < 0.99:
					status = "BAD!!!"
					result = "Fail"
					report = report + "\thas a correlation of " + correlations[ii] + ", below the accepted threshold of 0.99\n"

				RE = re.compile( '^Velocity' )
				RESULT = RE.match( keys[ii] )

				if RESULT != None:
					if cvgRates[ii] < 1.6:
						status = "BAD!!!"
						result = "Fail"
						report = report + "\thas a convergence rate of " + cvgRates[ii+1] + ", theoretically velocity fields should converge at a rate ~ 2\n"

				RE = re.compile( 'Pressure' )
				RESULT = RE.match( keys[ii] )

				if RESULT != None:
					if cvgRates[ii] < 0.9:
						status = "BAD!!!"
						result = "Fail"
						report = report + "\thas a convergence rate of " + cvgRates[ii+1] + ", theoretically pressure fields should converge at a rate ~ 1\n" 
				RE = re.compile( 'StrainRate' )
				RESULT = RE.match( keys[ii] )

				if RESULT != None:
					if cvgRates[ii] < 0.85:
						status = "BAD!!!"
						result = "Fail"
						report = report + "\thas a convergence rate of " + cvgRates[ii+1] + ", heoretically pressure fields should converge at a rate ~ 1\n"

				sys.stdout.write( status + "\n" + report )
		finally:
			INPUT.close()
	except IOError:
		pass
	
	# move cvg dat to the log directory
	command = "mv " + datFile_param + " log/"
	os.system( command )

	sys.stdout.write( "\nResult = " + result + "\n" )

runTests()
testConvergence( generateConvergence() )


