#! /usr/bin/env python

import os
import sys
import subprocess

def runTests():

    commands = ['./runAndTestConvergence.pl testVelicSolA.xml -optionsFile OFile2D.dat', \
                './runAndTestConvergence.pl testVelicSolB.xml -optionsFile OFile2D.dat', \
                './runAndTestConvergence.pl testVelicSolCx.xml -optionsFile OFile2D.dat', \
                './runAndTestConvergence.pl testVelicSolKz.xml -optionsFile OFile2D.dat', \
                './runAndTestConvergence.pl testVelicSolS.xml -optionsFile OFile2D.dat', \
                './runAndTestConvergence.pl testDepthDependentViscosity3D_Exponential.xml -optionsFile OFile3D.dat', \
                './sprTestAll.py', \
                './repTestAll.py' ]
    failed_commands = []
    passed = 0
    failed = 0

    for command in commands:
        try:
            retcode = subprocess.call( command, shell=True )
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
          "[SYS] Underworld Convergence Tests:\n" + \
          "[SYS]      Total Passes: (" + str(passed) + "/" + str(len( commands )) + ")\n" \

    if( len(failed_commands) > 0 ):
        message += "[SYS]      Failed Commands:\n"
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
