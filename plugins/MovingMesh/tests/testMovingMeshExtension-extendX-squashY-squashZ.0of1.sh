#!/bin/sh

TEST_SCRIPT=./VMake/executableTester.sh
until test -r ${TEST_SCRIPT} ; do
        TEST_SCRIPT=../${TEST_SCRIPT}
done
. ${TEST_SCRIPT}

runAndHandleSystemTestStdLocations "Underworld testExtension.xml velocityICsAndBCs.extensionX-squashY-squashZ.xml --remeshAccordingToJAxis=True --remeshAccordingToKAxis=True" "$0" "$@"
