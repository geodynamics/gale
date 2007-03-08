#Finds the Absolute path to the Project Root directory
SHELL := /bin/bash
PROJ_ROOT := $(shell until test -r ./Makefile.system ; do cd .. ; done ; echo `pwd`)
include ${PROJ_ROOT}/Makefile.system

# External Libraries and Headers
packages = PICELLERATOR STGFEM STGERMAIN MPI XML 

# Note don't put the .c extension. just the base name of the plugin
plugins := $(wildcard *.c)

include ${PROJ_ROOT}/Makefile.vmake

