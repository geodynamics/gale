#!/bin/sh

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
MaxTimeStepsDefault=3
SaveCheckpointStep=2
StepsLeft=1

# test if xml with this name exists
if test -r "${partTestname}.xml"; then
	# Run MPI's to generate results to test against.
	printf "$testname: doing pre-test setup:\n"
	
	printf "\tRunning for $MaxTimeStepsDefault timesteps to generate expected result\n"
	rm -rf ./output/"${partTestname}${procToWatch}of${nproc}GeneratedTestResult"
	RunMPICommand ${testname} Underworld "./${partTestname}.xml" --interactive=False --maxTimeSteps=$MaxTimeStepsDefault --checkpointEvery=1 --outputPath=./output/"${partTestname}${procToWatch}of${nproc}GeneratedTestResult" > ./log/"${testname}.generateExpectedResult.out" 2> ./log/"${testname}.generateExpectedResult.error"
	
	# copy checkpoint data to location to reload from for runAndHandleSystemTestStdLocations.
	printf "\tCopying results for $SaveCheckpointStep timesteps to reload location.\n"

	cp -rf ./output/"${partTestname}${procToWatch}of${nproc}GeneratedTestResult/"*"$SaveCheckpointStep".dat ./output/ 	
	# Do checkpointing test
	printf "Doing actual test: restarting from timestep $SaveCheckpointStep, checking if result by timestep $MaxTimeStepsDefault == previously generated one\n"
	runAndHandleSystemTestStdLocations "Underworld ./${partTestname}.xml --interactive=False --dumpEvery=1 --maxTimeSteps=$StepsLeft --restartTimestep=$SaveCheckpointStep --checkpointEvery=1 --plugins[]=StgFEM_CompareFeVariableAgainstReferenceSolution --StgFEM_CompareFeVariableAgainstReferenceSolution.referencePath=./output/${partTestname}${procToWatch}of${nproc}GeneratedTestResult --StgFEM_CompareFeVariableAgainstReferenceSolution.timeStepToCompare=$MaxTimeStepsDefault --StgFEM_CompareFeVariableAgainstReferenceSolution.referenceFeVariableSuffix=" "$0" "$@"

else echo "${partTestname}.xml doesn't exist, so can't test it"
fi
