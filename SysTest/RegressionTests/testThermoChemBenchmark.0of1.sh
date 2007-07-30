#!/bin/sh

TEST_SCRIPT=./VMake/executableTester.sh
until test -r ${TEST_SCRIPT} ; do
        TEST_SCRIPT=../${TEST_SCRIPT}
done
. ${TEST_SCRIPT}

setUpdateModeFromArg ${1}

testname="`basename $0 .sh`"
runAndHandleSystemTestStdLocations "Underworld testThermoChemBenchmark.xml --maxTimeSteps=5 --interactive=False --dumpEvery=5 --checkpointEvery=1" "$0" "$@"
