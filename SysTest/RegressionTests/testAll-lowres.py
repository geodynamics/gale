#! /usr/bin/env python

import os
import sys
import subprocess

def runTests():

    commands = ['./systest.pl Anisotropic.xml -optionsFile lowres_np1.dat -D lowres_expected', \
                  './systest.pl Arrhenius.xml -optionsFile lowres_np1.dat -D lowres_expected', \
                  './systest.pl ArrheniusPIC.xml -optionsFile lowres_np1.dat -D lowres_expected', \
                  './systest.pl CylinderRiseThermal.xml -optionsFile lowres_np1.dat -D lowres_expected', \
                  './systest.pl DepthDependentViscosity.xml -optionsFile lowres_np1.dat -D lowres_expected', \
                  './systest.pl Extension.xml -optionsFile lowres_np1.dat -D lowres_expected', \
                  './systest.pl FrankKamenetskii.xml -optionsFile lowres_np1.dat -D lowres_expected', \
                  './systest.pl InternalHeating.xml -optionsFile lowres_np1.dat -D lowres_expected', \
                  './systest.pl MobileLid.xml -optionsFile lowres_np1.dat -D lowres_expected', \
                  './systest.pl MultiThermalDiffusivity.xml -optionsFile lowres_np1.dat -D lowres_expected', \
                  './systest.pl MultiComponent.xml -optionsFile lowres_np1.dat -D lowres_expected', \
                  './systest.pl NonNewtonian.xml -optionsFile lowres_np1.dat -D lowres_expected', \
                  './systest.pl NonNewtonianPicard.xml -optionsFile lowres_np1.dat -D lowres_expected', \
                  './systest.pl RayleighTaylorBenchmark.xml -optionsFile lowres_np1.dat -D lowres_expected', \
                  './systest.pl ThermoChemBenchmark.xml -optionsFile lowres_np1.dat -D lowres_expected' ]
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
          "[SYS] Underworld Low-Res Integration Tests:\n" + \
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
