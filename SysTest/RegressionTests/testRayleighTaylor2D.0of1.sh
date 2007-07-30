#!/bin/sh

TEST_SCRIPT=./VMake/executableTester.sh
until test -r ${TEST_SCRIPT} ; do
        TEST_SCRIPT=../${TEST_SCRIPT}
done
. ${TEST_SCRIPT}

setUpdateModeFromArg ${1}

testname="`basename $0 .sh`"
runAndHandleSystemTestStdLocations "Underworld testRayleighTaylorBenchmark.xml --interactive=False --maxTimeSteps=5 --checkpointEvery=1 --plugins[]=StgFEM_CompareFeVariableAgainstReferenceSolution --StgFEM_CompareFeVariableAgainstReferenceSolution.referencePath=./expected/testRayleighTaylor2D.0of1" "$0" "$@"
