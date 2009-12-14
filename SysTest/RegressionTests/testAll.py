#! /usr/bin/env python

import os

def runTests():

    commands = ['./analyticTest.pl CosineHillRotateBC.xml -optionsFile np-1.dat', \
                './analyticTest.pl CosineHillRotateBC.xml -optionsFile np-2.dat', \
                './analyticTest.pl CosineHillRotateBC.xml -optionsFile np-4.dat', \
                './analyticTest.pl CosineHillRotateBC-DualMesh.xml -optionsFile np-1.dat', \
                './analyticTest.pl CosineHillRotateBC-DualMesh.xml -optionsFile np-2.dat', \
                './analyticTest.pl CosineHillRotateBC-DualMesh.xml -optionsFile np-4.dat', \
                './analyticTest.pl HomogeneousNaturalBCs.xml -optionsFile np-1.dat', \
                './analyticTest.pl HomogeneousNaturalBCs.xml -optionsFile np-2.dat', \
                './analyticTest.pl HomogeneousNaturalBCs.xml -optionsFile np-4.dat', \
                './analyticTest.pl HomogeneousNaturalBCs-DualMesh.xml -optionsFile np-1.dat', \
                './analyticTest.pl HomogeneousNaturalBCs-DualMesh.xml -optionsFile np-2.dat', \
                './analyticTest.pl HomogeneousNaturalBCs-DualMesh.xml -optionsFile np-4.dat', \
                './analyticTest.pl SteadyState1D-x.xml -optionsFile ss-0.5-np-1.dat', \
                './analyticTest.pl SteadyState1D-x.xml -optionsFile ss-0.5-np-2.dat', \
                './analyticTest.pl SteadyState1D-x.xml -optionsFile ss-0.5-np-4.dat', \
                './analyticTest.pl SteadyState1D-x.xml -optionsFile np-1.dat', \
                './analyticTest.pl SteadyState1D-x.xml -optionsFile np-2.dat', \
                './analyticTest.pl SteadyState1D-x.xml -optionsFile np-4.dat', \
                './analyticTest.pl SteadyState1D-y.xml -optionsFile np-1.dat', \
                './analyticTest.pl SteadyState1D-y.xml -optionsFile np-2.dat', \
                './analyticTest.pl SteadyState1D-y.xml -optionsFile np-4.dat', \
                './analyticTest.pl AnalyticSimpleShear.xml -optionsFile np-1.dat', \
                './analyticTest.pl AnalyticSimpleShear.xml -optionsFile np-2.dat', \
                './analyticTest.pl AnalyticSimpleShear.xml -optionsFile np-4.dat', \
                './analyticTest.pl AnalyticSinusoid.xml -optionsFile as-np-1.dat', \
                './analyticTest.pl AnalyticSinusoid.xml -optionsFile as-np-2.dat', \
                './analyticTest.pl AnalyticSinusoid.xml -optionsFile as-np-4.dat', \
                './analyticTest.pl TempDiffusion.xml -optionsFile np-1.dat', \
                './analyticTest.pl TempDiffusion.xml -optionsFile np-2.dat', \
                './analyticTest.pl TempDiffusion.xml -optionsFile np-4.dat', \
                './analyticTest.pl TempDiffusion.xml -optionsFile np-8.dat']

    for command in commands:
        os.system( command )

runTests()
