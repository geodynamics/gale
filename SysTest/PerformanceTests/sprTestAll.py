#! /usr/bin/env python

import os

def runTests():

    commands = ['./runAndTestConvergence.pl testVelicSolA.xml -optionsFile OFileSPR.dat', \
                './runAndTestConvergence.pl testVelicSolB.xml -optionsFile OFileSPR.dat', \
                './runAndTestConvergence.pl testVelicSolCx.xml -optionsFile OFileSPR.dat', \
                './runAndTestConvergence.pl testVelicSolKz.xml -optionsFile OFileSPR.dat', \
                './runAndTestConvergence.pl testVelicSolS.xml -optionsFile OFileSPR.dat' ]

    for command in commands:
        os.system( command )


runTests()
