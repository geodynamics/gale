#!/bin/sh

# Input parameters
expectedCleaned=$1
savedResultCleaned=$2

# NB - Should get value from Makefile.system
DIFF=diff

# Run 'diff' on the files with the opengl statements in them
diffs=$savedResultCleaned.diffs.tmp
${DIFF} ${expectedCleaned} ${savedResultCleaned} > ${diffs}

# Run a perl script to weed out all the small numerical differences
./CleanOpenGLDiffFile.pl ${diffs} > ${diffs}.cleaned
mv ${diffs}.cleaned ${diffs}

# Test if $diffs is not an empty file - if it isn't then this test failed
if test -s ${diffs} ; then
	cat ${diffs}
	rm ${diffs}
	exit 1
fi

# if $diffs is an empty file then this test passed
rm $diffs
exit 0
