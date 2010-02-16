#! /usr/bin/env python

import os
import sys
import subprocess

def runTests():

    commands = ['./systest.pl Anisotropic.xml -n 10', \
                './systest.pl Arrhenius.xml -n 10', \
                './systest.pl Arrhenius.xml -optionsFile np-2.dat -n 10', \
                './systest.pl ArrheniusPIC.xml -n 10', \
                './systest.pl CylinderRiseThermal.xml -n 10 ', \
                './systest.pl DepthDependentViscosity.xml -n 10', \
                './systest.pl Extension.xml -n 10', \
                './systest.pl FrankKamenetskii.xml -n 10', \
                './systest.pl InternalHeating.xml -n 10', \
                './systest.pl MobileLid.xml -n 10', \
                './systest.pl MultiThermalDiffusivity.xml -n 10', \
                './systest.pl MultiComponent.xml -n 10', \
                './systest.pl NonNewtonian.xml -n 10', \
                './systest.pl NonNewtonianPicard.xml -n 10', \
                './systest.pl RayleighTaylorBenchmark.xml -n 10', \
                './systest.pl RayleighTaylorBenchmark.xml -optionsFile np-2.dat -n 10', \
                './systest.pl ThermoChemBenchmark.xml -n 10', \
                './analyticTest.pl NonNewtonianShear.xml -optionsFile np-1.dat', \
                './analyticTest.pl NonNewtonianShear.xml -optionsFile np-2.dat', \
                './analyticTest.pl NonNewtonianShear.xml -optionsFile np-4.dat', \
                './analyticTest.pl NonNewtonianShear.xml -optionsFile np-8.dat', \
                './analyticTest.pl Trubitsyn2006Isoviscous.xml -optionsFile np-1.dat', \
                './analyticTest.pl Trubitsyn2006Isoviscous.xml -optionsFile np-2.dat', \
                './analyticTest.pl Trubitsyn2006Isoviscous.xml -optionsFile np-4.dat', \
                './analyticTest.pl Trubitsyn2006Isoviscous.xml -optionsFile np-8.dat', \
                './restartTest.pl Arrhenius.xml', \
                './restartTest.pl DepthDependentViscosity.xml', \
                './restartTest.pl Extension.xml', \
                './restartTest.pl NonNewtonian.xml', \
                './restartTest.pl RayleighTaylorBenchmark.xml' ]
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
          "[SYS] Underworld Normal-Res Integration Tests:\n" + \
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
