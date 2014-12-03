#!/bin/sh
#Find testing shell scripts
TEST_SCRIPT=./VMake/executableTester.sh
until test -r ${TEST_SCRIPT} ; do
        TEST_SCRIPT=../${TEST_SCRIPT}
done
. ${TEST_SCRIPT}
. `dirname ${TEST_SCRIPT}`/build-functions.sh

setUpdateModeFromArg ${1}

# Extract out names
testname="`basename $0 .sh`"
nproc=`echo ${testname} | cut -d . -f 2 | cut -c 4`
procToWatch=`echo ${testname} | cut -d . -f 2 | cut -f 1 -d 'o'`
partTestname=`echo $testname | cut -f 1 -d '-'`

# test if xml with this name exists
if test -r "${partTestname}.xml"; then
	# Run MPI's to generate results to test against.
	printf "$testname: doing pre-test setup:\n"
	
	printf "\tRunning for 4 timesteps to generate expected result\n"
	rm -rf ./output/"${partTestname}${procToWatch}of${nproc}GeneratedTestResult"
	RunMPICommand ${testname} Underworld "./${partTestname}.xml" --interactive=False --maxTimeSteps=4 --checkpointEvery=1 --elementResI=24 --elementResJ=24 --elementResK=3 --outputPath=./output/"${partTestname}${procToWatch}of${nproc}GeneratedTestResult" > ./log/"${testname}.generateExpectedResult.out" 2> ./log/"${testname}.generateExpectedResult.error"
	
	# Run MPI's to generate first few stepa to checkpoint.
	printf "\tRunning for 2 timesteps to generate checkpoints to reload from\n"
	RunMPICommand ${testname} Underworld "./${partTestname}.xml" --interactive=False --maxTimeSteps=2 --checkpointEvery=1 --elementResI=24 --elementResJ=24 --elementResK=3  > ./log/"${testname}.generateCheckpoint.out" --outputPath=./output/ 2> ./log/"${testname}.generateCheckpoint.error"
	
	# Do checkpointing test
	printf "Doing actual test: restarting from timestep 2, checking if result by timestep 4 == previously generated one\n"
	runAndHandleSystemTestStdLocations "Underworld ./${partTestname}.xml --interactive=False --dumpEvery=1 --maxTimeSteps=2 --restartTimestep=2 --checkpointEvery=1 --elementResI=24 --elementResJ=24 --elementResK=3 --plugins[]=StgFEM_CompareFeVariableAgainstReferenceSolution --StgFEM_CompareFeVariableAgainstReferenceSolution.referencePath=./output/${partTestname}${procToWatch}of${nproc}GeneratedTestResult --StgFEM_CompareFeVariableAgainstReferenceSolution.timeStepToCompare=4 --StgFEM_CompareFeVariableAgainstReferenceSolution.referenceFeVariableSuffix=" "$0" "$@"

else echo "${partTestname}.xml doesn't exist, so can't test it"
fi
