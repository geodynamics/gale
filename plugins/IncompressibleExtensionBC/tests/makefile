
#Finds the Absolute path to the Project Root directory
SHELL := /bin/bash
PROJ_ROOT := $(shell until test -r ./Makefile.system ; do cd .. ; done ; echo `pwd`)
include ${PROJ_ROOT}/Makefile.system

# Subdirectories
subdirs := IncompressibleExtensionAnalytic

checks = testIncompressibleExtensionBC.0of1.sh \
	testIncompressibleExtensionBC3D.0of1.sh

long_checks = testIncompressibleExtensionBC.0of4.sh \
	testIncompressibleExtensionBC3D.0of4.sh

include ${PROJ_ROOT}/Makefile.vmake
