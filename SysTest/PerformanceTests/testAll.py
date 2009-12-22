#! /usr/bin/env python

import os

def runTests():

    commands = ['./runAndTestConvergence.pl testVelicSolA.xml -optionsFile OFile2D.dat', \
                './runAndTestConvergence.pl testVelicSolB.xml -optionsFile OFile2D.dat', \
                './runAndTestConvergence.pl testVelicSolCx.xml -optionsFile OFile2D.dat', \
                './runAndTestConvergence.pl testVelicSolKz.xml -optionsFile OFile2D.dat', \
                './runAndTestConvergence.pl testVelicSolS.xml -optionsFile OFile2D.dat', \
                './runAndTestConvergence.pl testDepthDependentViscosity3D_Exponential.xml -optionsFile OFile3D.dat', \
                './sprTestAll.py', \
                './repTestAll.py' ]

    for command in commands:
        os.system( command )


runTests()
