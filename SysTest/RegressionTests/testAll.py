#! /usr/bin/env python

import os

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

    for command in commands:
        os.system( command )

runTests()







