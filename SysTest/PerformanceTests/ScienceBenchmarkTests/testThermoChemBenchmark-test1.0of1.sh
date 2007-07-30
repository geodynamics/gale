#!/bin/sh

#########################################################################
# FILE PURPOSE:
#
# This test is designed to be a framework for a scientific benchmark test.
# It will load checkpointed values from an expected directory specific
# for the test and machine type it is testing on,
# Run the xml name given in the title of the test ( minus -test1.0of1.sh)
# Use the python script in the plugins directory to read the FrequentOutput.dat 
# file for the specific valueName given, and compare the results from a
# given time (with tolerance) to an expected result for this value.
# It will ten write out a ${valueName}-result.dat file that will
# say whether the test passed or failed.
# A shell script function will compare this result.dat file against
# an expected file in the baseExpectedPath, and if they are the same
# the test passes.
###########################################################################

###########################################################################
# TEST CREATOR INSTRUCTIONS:
#
# To create your own tests:
# * Create an xml file with the name: ${NAMEOFXMLFILE}.xml if not using
#   an xml file already in the directory
# * Copy this shell script file to another name 
# 	The name MUST have the following format:
#         ${NAMEOFXMLFILE}-${TESTPREFIX}.${PROCTOWATCH}of${NPROC}.sh 
# * Create a directory in ./expected called test${NAMEOFXMLFILE}
# * Create a directory for each machine-type the test will be run on.
#   The machine-types for each cluster can be found in the appropriate 
#	Makefile.system file under MACHINES=
#	eg. local machines are often MACGINES=i686
# * Put in checkpointed data for for each $MACHINES sub-dir for the 
#	checkpointed step to run from. 
# * Set CheckpointValue to load from in this file.
# * Set MaxTimeStepsDefault in this file so program knows how many time steps to run
#	from loaded checkpoint timestep.
# * Go to middle of this file and change scientific values:
#			valueName
#			expectedValue
#			toleranceValue
#			expectedTime
#			toleranceTime
#	Read explanation to know what format these values will need to be in.
# * Change the paper reference so that the user will know what paper
#	the test is benchmarked against.
# * Create an expected file in ./expected/test${NAMEOFXMLFILE} with the 
#	name format: ${TESTPREFIX}-${valueName}-result.out.expected
#	This file will need to be patched before it will be accurate.
# * To (dodgily) patch a file: run the test, and copy the test result from
#	./output/test${NAMEOFXMLFILE}/$MACHINES/${TESTPREFIX}-${valueName}-result.out
#	to the expected file.
#	Not sure if the standard patch method will work at present.
# * To add test to regression system:
#		go to ScienceBenchMarkTests directory ( if that is where your test is )
#		edit the makefile, under 	scibenchmark_checks =  \
#		add in the line (with your own names for each of the ${} entries):
#		test${NAMEOFXMLFILE}-checkpoint-${TESTPREFIX}.${PROCTOWATCH}of${NPROC}.sh \
#		save the make file.
# * To test your shell script file, either run:
#		- make scibenchmark_check 
#		- sh test${NAMEOFXMLFILE}-checkpoint-${TESTPREFIX}.${PROCTOWATCH}of${NPROC}.sh
#		- ./test${NAMEOFXMLFILE}-checkpoint-${TESTPREFIX}.${PROCTOWATCH}of${NPROC}.sh
#	MAKE SURE THE SHELL SCRIPT IS EXECUTABLE (use chmod to change it if it isn't)
# * Only when the test is working and passing should you commit your test to the repository.
###################################################################################

#################
# TEST:  



#### CHANGE THESE VALUES ####
# This section sets the values to run from and to.
# CheckpointValue is the timestep that the program will
# be loaded from
CheckpointValue=3
# MaxTimeStepsDefault is the number of timesteps that 
# will be run from the checkpointed data
MaxTimeStepsDefault=2
#############################

# This recursion assumes program is a StGermain-style project with
# a VMake directory.
TEST_SCRIPT=./VMake/executableTester.sh
until test -r ${TEST_SCRIPT} ; do
        TEST_SCRIPT=../${TEST_SCRIPT}
done
. ${TEST_SCRIPT}
. `dirname ${TEST_SCRIPT}`/build-functions.sh

setUpdateModeFromArg ${1}

# Grab system information from Makefile.system
getValueFromMakefile_System MACHINE
getValueFromMakefile_System BIN_DIR
getValueFromMakefile_System MPI_NPROC
getValueFromMakefile_System MPI_RUN
getValueFromMakefile_System MPI_SGIIMPLEMENT
getValueFromMakefile_System MPI_MACHINES


# This section extracts out the testnames, procToWatch and nproc
# from the shell script name
testname="`basename $0 .sh`"
nproc=`echo ${testname} | cut -d . -f 2 | cut -c 4`
procToWatch=`echo ${testname} | cut -d . -f 2 | cut -f 1 -d 'o'`
partTestname=`echo $testname | cut -f 1 -d '-'`

# testExtension is used to differentiate between this test 
# and one for a different checkpointed data test using the same xml.	

testExtension=`echo ${testname} | cut -f 2 -d '-' | cut -d . -f 1`

# This section sets the where the outputfiles are to be printed
# and where the expected files and checkpointed data are to be found.
outputPath=./output/${testname}/$MACHINE/
expectedPath=./expected/${testname}/$MACHINE/
baseOutputPath=./output/${testname}
baseExpectedPath=./expected/${testname}

# Test if expectedPath and baseExpectedPath exist.
if test -d ${baseExpectedPath} ; then
	if test -d ${expectedPath} ; then
		# test if xml with this name exists
		if test -r "${partTestname}.xml"; then
			# copy checkpointed data from correct ($MACHINES) expected directory
			printf "Testing Scientific Benchmark: ${testname}:\n"
			printf "\tBenchmark results taken from reference:\n"
########### PUT IN REFERENCE TO PAPER HERE ##################
			# Here is where you put your reference to the scientific paper
			# This test is comparing against.
			printf "\t\tPUT YOUR REFERENCE TO PAPER HERE \n"
#############################################################			
			if test -d ${outputPath}; then
				printf "\tRemoving old data from $outputPath\n"
				rm -rf ${outputPath}/*
			else 
				if ! test -d ${baseOutputPath}; then		
					printf "\tCreating base test directory $baseOutputPath\n"
					mkdir ${baseOutputPath}
				fi
				printf "\tCreating output directory $outputPath\n"
				mkdir ${outputPath}		
			fi
			
			printf "\tCopying info for timestep = $CheckpointValue from expected directory $expectedPath to outputPath=$outputPath\n"
			cp ${expectedPath}/* ${outputPath}
			printf "\tRunning from timestep $CheckpointValue for $MaxTimeStepsDefault timesteps to generate expected result\n"
			# Run MPI to generate FrequentOutput.dat file
			${MPI_RUN} ${MPI_SGIIMPLEMENT} ${MPI_MACHINES} ${MPI_NPROC} ${nproc} ${BIN_DIR}/Underworld ./${partTestname}.xml --interactive=False --dumpEvery=1 --maxTimeSteps=$MaxTimeStepsDefault --restartTimestep=$CheckpointValue --outputPath=${outputPath} > ./log/"${testname}.out" 2> ./log/"${testname}.error"
			
			# Test to see if Frequentoutput.dat is of required level of accuracy
			# using python script in plugins directory:
			# From python help file:
			# usage: readFrequentOutput.py [options]
			#options:
			#-h, --help            show this help message and exit
			#-o OUTPUTPATH, --outputPath=OUTPUTPATH
			#                    output Directory path
			#-n VALUENAME, --valueName=VALUENAME
			#                    Name of value in FrequentOutput.dat file header
			#-V EXPECTEDVALUE, --expectedValue=EXPECTEDVALUE
			#                    expected numerical Value at given time
			#-T EXPECTEDTIME, --expectedTime=EXPECTEDTIME
			#                    expected Time value
			#-t TOLERANCETIME, --toleranceTime=TOLERANCETIME
			#                    tolerance on time value
			#-v TOLERANCEVALUE, --toleranceValue=TOLERANCEVALUE
			#                    tolerance on expected value at expected time
		
			

#### CHANGE THESE VALUES ####
			# Set value to be tested
			# These values will need to change to match scientific benchmark values
			# for each test.
			valueName=Vrms
			expectedValue=0.0636
			toleranceValue=0.01
			expectedTime=0.00073
			toleranceTime=0.00001
################################
			
			# Create the results output filename from the known
			# inputs to the python script.
			resultsFilename=${testExtension}-${valueName}-result.out
			
			# Remove old result files
			rm -rf ${outputPath}/${resultsFilename}
		
			printf "\tChecking to see if at time ${expectedTime} +/- ${toleranceTime}, value of ${valueName} is ${expectedValue} +/- ${toleranceValue} as expected from (ref)\n"
			python plugins/readFrequentOutput.py --outputPath=${outputPath} --valueName=${valueName} --expectedValue=${expectedValue} --toleranceValue=${toleranceValue} --expectedTime=${expectedTime} --toleranceTime=${toleranceTime}	
			
			mv ${outputPath}/${valueName}-result.out ${outputPath}/${resultsFilename}
			
			# Now run test to see if the results match the expected files
################################################################################
#					DO NOT TOUCH THE CODE BELOW. 
#
#					MAKE SURE IT IS IN EACH TEST FILE
################################################################################
			
			# TODO This is still a hack: see executableTester.sh 
			RunScienceTestCheckStdLocations "${testname}" "${resultsFilename}" 
		
		else echo "${partTestname}.xml doesn't exist, so can't test it"
		
		fi
		
		# Lastly, remove checkpointed data from outputPath
		rm -rf ${outputPath}/*.dat
	else
		printf "${expectedPath} does not exist. Cannot run test, ${testname} for machine = $MACHINE \n"
	fi
else
	printf "${baseExpectedPath} does not exist. Cannot run any tests for ${testname} \n"
fi
