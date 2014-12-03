#!/bin/sh

TEST_SCRIPT=./VMake/executableTester.sh
until test -r ${TEST_SCRIPT} ; do
        TEST_SCRIPT=../${TEST_SCRIPT}
done
. ${TEST_SCRIPT}

setUpdateModeFromArg ${1}

testname="`basename $0 .sh`"
runAndHandleSystemTestStdLocations "Underworld testBlankbachCase1a.xml "$0" "$@"
