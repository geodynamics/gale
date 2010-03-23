#! /usr/bin/env python

import os
import sys
import subprocess

def runTests():

    commands = ['./analyticTest.pl CosineHillRotateBC.xml -optionsFile np-1.dat', \
                './analyticTest.pl CosineHillRotateBC.xml -optionsFile np-2.dat', \
                './analyticTest.pl CosineHillRotateBC.xml -optionsFile np-4.dat', \
                './analyticTest.pl CosineHillRotateBC-DualMesh.xml -optionsFile np-1.dat', \
                './analyticTest.pl CosineHillRotateBC-DualMesh.xml -optionsFile np-2.dat', \
                './analyticTest.pl CosineHillRotateBC-DualMesh.xml -optionsFile np-4.dat', \
                './analyticTest.pl HomogeneousNaturalBCs.xml -optionsFile np-1.dat', \
                './analyticTest.pl HomogeneousNaturalBCs.xml -optionsFile np-2.dat', \
                './analyticTest.pl HomogeneousNaturalBCs.xml -optionsFile np-4.dat', \
                './analyticTest.pl HomogeneousNaturalBCs-DualMesh.xml -optionsFile np-1.dat', \
                './analyticTest.pl HomogeneousNaturalBCs-DualMesh.xml -optionsFile np-2.dat', \
                './analyticTest.pl HomogeneousNaturalBCs-DualMesh.xml -optionsFile np-4.dat', \
                './analyticTest.pl SteadyState1D-x.xml -optionsFile ss-0.5-np-1.dat', \
                './analyticTest.pl SteadyState1D-x.xml -optionsFile ss-0.5-np-2.dat', \
                './analyticTest.pl SteadyState1D-x.xml -optionsFile ss-0.5-np-4.dat', \
                './analyticTest.pl SteadyState1D-x.xml -optionsFile np-1.dat', \
                './analyticTest.pl SteadyState1D-x.xml -optionsFile np-2.dat', \
                './analyticTest.pl SteadyState1D-x.xml -optionsFile np-4.dat', \
                './analyticTest.pl SteadyState1D-y.xml -optionsFile np-1.dat', \
                './analyticTest.pl SteadyState1D-y.xml -optionsFile np-2.dat', \
                './analyticTest.pl SteadyState1D-y.xml -optionsFile np-4.dat', \
                './analyticTest.pl AnalyticSimpleShear.xml -optionsFile np-1.dat', \
                './analyticTest.pl AnalyticSimpleShear.xml -optionsFile np-2.dat', \
                './analyticTest.pl AnalyticSimpleShear.xml -optionsFile np-4.dat', \
                './analyticTest.pl AnalyticSinusoid.xml -optionsFile as-np-1.dat', \
                './analyticTest.pl AnalyticSinusoid.xml -optionsFile as-np-2.dat', \
                './analyticTest.pl AnalyticSinusoid.xml -optionsFile as-np-4.dat', \
                './analyticTest.pl TempDiffusion.xml -optionsFile np-1.dat', \
                './analyticTest.pl TempDiffusion.xml -optionsFile np-2.dat', \
                './analyticTest.pl TempDiffusion.xml -optionsFile np-4.dat', \
                './analyticTest.pl TempDiffusion.xml -optionsFile np-8.dat', \
                './checkpointTest.pl Multigrid.xml']
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
          "[SYS] StgFEM Normal-Res Integration Tests:\n" + \
          "[SYS]      Total Passes: (" + str(passed) + "/" + str(len( commands )) + ")\n"

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
