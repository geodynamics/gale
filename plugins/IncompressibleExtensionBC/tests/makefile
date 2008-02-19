
#Finds the Absolute path to the Project Root directory
SHELL := /bin/bash
PROJ_ROOT := $(shell until test -r ./Makefile.system ; do cd .. ; done ; echo `pwd`)
include ${PROJ_ROOT}/Makefile.system

# Subdirectories
subdirs := IncompressibleExtensionAnalytic

checks = $(wildcard *.sh)

include ${PROJ_ROOT}/Makefile.vmake
