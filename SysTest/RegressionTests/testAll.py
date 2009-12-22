#! /usr/bin/env python

import os

def runTests():

    commands = ['./checkpointTest.pl Anisotropic.xml', \
                './checkpointTest.pl Arrhenius.xml', \
                './checkpointTest.pl Arrhenius.xml -optionsFile np-2.dat', \
                './checkpointTest.pl ArrheniusPIC.xml', \
                './checkpointTest.pl CylinderRiseThermal.xml ', \
                './checkpointTest.pl DepthDependentViscosity.xml', \
                './checkpointTest.pl Extension.xml', \
                './checkpointTest.pl FrankKamenetskii.xml', \
                './checkpointTest.pl InternalHeating.xml', \
                './checkpointTest.pl MobileLid.xml', \
                './checkpointTest.pl MultiThermalDiffusivity.xml', \
                './checkpointTest.pl MultiComponent.xml', \
                './checkpointTest.pl NonNewtonian.xml', \
                './checkpointTest.pl NonNewtonianPicard.xml', \
                './checkpointTest.pl RayleighTaylorBenchmark.xml', \
                './checkpointTest.pl RayleighTaylorBenchmark.xml -optionsFile np-2.dat', \
                './checkpointTest.pl ThermoChemBenchmark.xml', \
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







