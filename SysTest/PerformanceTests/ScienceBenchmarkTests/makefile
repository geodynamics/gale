
#Finds the Absolute path to the Project Root directory
SHELL := /bin/bash
PROJ_ROOT := $(shell until test -r ./Makefile.system ; do cd .. ; done ; echo `pwd`)
include ${PROJ_ROOT}/Makefile.system


subdirs := plugins

scibenchmark_checks = testThermoChemBenchmark-test1.0of1.sh 

# These tests only work on 32 bit machines since that's where we generated the 
#	checkpoint "benchmark" files, and as yet the Swarm checkpoint files are
#	binary & not portable between architectures
	
ifeq (${SYSTEM_64},0)
	scibenchmark_checks += 

endif


include ${PROJ_ROOT}/Makefile.vmake
