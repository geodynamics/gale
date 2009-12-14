#! /usr/bin/env python

import os

def runTests():

    commands = ['./runAndTestConvergence.pl AnalyticSinusoid.xml -optionsFile as-np-1.dat', \
                './runAndTestConvergence.pl AnalyticSinusoid.xml -optionsFile as-np-2.dat', \
                './runAndTestConvergence.pl AnalyticSinusoid.xml -optionsFile as-np-4.dat']

    for command in commands:
        os.system( command )

runTests()
