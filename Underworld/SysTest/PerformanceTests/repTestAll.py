#! /usr/bin/env python

import os

def runTests():

    commands = ['./runAndTestConvergence.pl testVelicSolA.xml -optionsFile OFileREP.dat', \
                './runAndTestConvergence.pl testVelicSolB.xml -optionsFile OFileREP.dat', \
                './runAndTestConvergence.pl testVelicSolCx.xml -optionsFile OFileREP.dat', \
                './runAndTestConvergence.pl testVelicSolKz.xml -optionsFile OFileREP.dat', \
                './runAndTestConvergence.pl testVelicSolS.xml -optionsFile OFileREP.dat']

    for command in commands:
        os.system( command )


runTests()
