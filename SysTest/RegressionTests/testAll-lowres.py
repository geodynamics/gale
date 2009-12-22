#! /usr/bin/env python

import os

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

    for command in commands:
        os.system( command )

runTests()







