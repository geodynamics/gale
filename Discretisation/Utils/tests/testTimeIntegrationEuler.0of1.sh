#!/bin/sh

TEST_SCRIPT=./VMake/executableTester.sh
until test -r ${TEST_SCRIPT} ; do
        TEST_SCRIPT=../${TEST_SCRIPT}
done
. ${TEST_SCRIPT}

runAndHandleSystemTest "testTimeIntegration --DerivName0=Constant --DerivName1=Constant2 --order=1" "$0" "$@"
