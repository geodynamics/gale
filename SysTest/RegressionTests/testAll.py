#! /usr/bin/env python

import os
import sys
import subprocess

def runTests():
    
    commands = ['./analyticTest.pl AnalyticShearXZ.xml -optionsFile np-1.dat', \
                './analyticTest.pl AnalyticShearXZ.xml -optionsFile np-2.dat', \
                './analyticTest.pl AnalyticShearXZ.xml -optionsFile np-4.dat', \
                './analyticTest.pl AnalyticShearXZ.xml -optionsFile np-8.dat']

    failed_commands = [] 
    passed = 0
    failed = 0
 
    for command in commands:
        try:
            retcode = subprocess.call( command+' -serial' , shell=True )
            if retcode == 0:
                 passed += 1
            else:
                 failed += 1
                 failed_commands.append( command )
        except OSError, e:
            print >>sys.stderr, "Execution Failed:", e

    filename = "../../../summary.dat"

    if os.path.exists( filename ):
        FILE = open( filename, "a" )
    else:
        FILE = open( filename, "w" )

    message = ''
    message += "--------------------------------------------------------\n" + \
          "[SYS] PICellerator Normal-Res Integration Tests:\n" + \
          "[SYS]      Total Passes: (" + str(passed) + "/" + str(len( commands )) + ")\n" \
          "[SYS]      Failed Commands:\n"
    for command in failed_commands:
        message += "[SYS]            " + command + "\n"
    message += "--------------------------------------------------------\n"
    FILE.write( message )
    print message
    FILE.close()

    if failed > 0:
        sys.exit(1)
    else:
        sys.exit(0)

runTests()
