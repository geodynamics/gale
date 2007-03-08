
#Finds the Absolute path to the Project Root directory
SHELL := /bin/bash
PROJ_ROOT := $(shell until test -r ./Makefile.system ; do cd .. ; done ; echo `pwd`)
include ${PROJ_ROOT}/Makefile.system

# The checkpointing tests can be run on any system, since all the checkpoints are generated
#	dynamically
# Moved all the serial checkpointing to the front as these tests are a lot quicker
# than the parallel tests. They should be left in this order so that an overall system test
# can be done first before embarking on parallel checkpointing tests --Kathleen 11/10/2006
checks = testAnisotropic-checkpoint.0of1.sh \
	testArrhenius-checkpoint.0of1.sh \
	testExtension-checkpoint.0of1.sh \
	testExtensionWithPop.0of1.sh \
	testExtensionFMM-checkpoint.0of1.sh \
	testFrankKamenetskii-checkpoint.0of1.sh \
	testInternalHeating-checkpoint.0of1.sh \
	testMobileLid-checkpoint.0of1.sh \
	testNonNewtonian-checkpoint.0of1.sh \
	testMultiComponent-checkpoint.0of1.sh \
	testRayleighTaylor2D-checkpoint.0of1.sh \
	testSlabSubduction-checkpoint.0of1.sh \
	testThermoChemBenchmark-checkpoint.0of1.sh \
	testExtensionWithPopControl-checkpoint.0of1.sh
	
long_checks = testAnisotropic-checkpoint.0of2.sh \
	testArrhenius-checkpoint.0of2.sh \
	testExtension-checkpoint.0of2.sh \
	testFrankKamenetskii-checkpoint.0of2.sh \
	testExtensionFMM-checkpoint.0of2.sh \
	testInternalHeating-checkpoint.0of2.sh \
	testMobileLid-checkpoint.0of2.sh \
	testMultiComponent-checkpoint.0of2.sh \
	testNonNewtonian-checkpoint.0of2.sh \
	testRayleighTaylor2D-checkpoint.0of2.sh \
	testSlabSubduction-checkpoint.0of2.sh \
	testThermoChemBenchmark-checkpoint.0of2.sh \
	testAnisotropic-checkpoint.0of4.sh \
	testArrhenius-checkpoint.0of4.sh \
	testFrankKamenetskii-checkpoint.0of4.sh \
	testInternalHeating-checkpoint.0of4.sh \
	testMobileLid-checkpoint.0of4.sh \
	testMultiComponent-checkpoint.0of4.sh \
	testNonNewtonian-checkpoint.0of4.sh \
	testSlabSubduction-checkpoint.0of4.sh \
	testThermoChemBenchmark-checkpoint.0of4.sh

# These tests only work on 32 bit machines since that's where we generated the 
#	checkpoint "benchmark" files, and as yet the Swarm checkpoint files are
#	binary & not portable between architectures
	
ifeq (${SYSTEM_64},0)
	checks += \
		testExtensionFMM.0of1.sh \
		testRayleighTaylor2D.0of1.sh \
		testThermoChemBenchmark.0of1.sh 
 
	long_checks += \
		testExtensionFMM.0of2.sh \
		testRayleighTaylor2D.0of2.sh \
		testThermoChemBenchmark.0of2.sh		
endif

#TODO Investigating why these tests are failing:
	# "particle outside box" at 30x15 - why?
	#testExtensionDP.0of1.sh \
	#testExtensionDP-checkpoint.0of1.sh \
	#testExtensionFMM3D-checkpoint.0of1.sh \
#subdirs := 

include ${PROJ_ROOT}/Makefile.vmake
