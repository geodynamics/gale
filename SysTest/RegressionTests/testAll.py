#! /usr/bin/env python

import os

def runTests():
    
    commands = ['./analyticTest.pl AnalyticShearXZ.xml -optionsFile np-1.dat', \
                './analyticTest.pl AnalyticShearXZ.xml -optionsFile np-2.dat', \
                './analyticTest.pl AnalyticShearXZ.xml -optionsFile np-4.dat', \
                './analyticTest.pl AnalyticShearXZ.xml -optionsFile np-8.dat']

    for command in commands:
        os.system( command )

runTests()
