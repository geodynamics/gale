#Finds the Absolute path to the Project Root directory
SHELL := /bin/bash
PROJ_ROOT := $(shell until test -r ./Makefile.system ; do cd .. ; done ; echo `pwd`)
include ${PROJ_ROOT}/Makefile.system

checks = testVelicSolA.0of1.sh \
	 testVelicSolB.0of1.sh \
	 testVelicSolCx.0of1.sh \
	 testVelicSolHA.0of1.sh \
	 testVelicSolKz.0of1.sh \
	 testVelicSolS.0of1.sh \

# Subdirectories
subdirs := ScienceBenchmarkTests

include ${PROJ_ROOT}/Makefile.vmake
