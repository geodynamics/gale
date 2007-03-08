#Finds the Absolute path to the Project Root directory
SHELL := /bin/bash
PROJ_ROOT := $(shell until test -r ./Makefile.system ; do cd .. ; done ; echo `pwd`)
include ${PROJ_ROOT}/Makefile.system

# Source Code and Header files
SRCS := $(wildcard *.c)
HDRS := $(wildcard *.h)

# External Libraries and Headers
EXTERNAL_INCLUDES = -DLIB_DIR=\"${LIB_DIR}\"
packages = PICELLERATOR STGFEM STGERMAIN PETSC MPI XML PYTHON

# Work out name to call library
lib = ${PROJECT}$(shell basename `pwd | sed s/src//g`)
libdynamic = yes
libstatic = yes

# Name of directory to place header files
includes := ${PROJECT}/$(shell basename `pwd | sed s/src//g`)

include ${PROJ_ROOT}/Makefile.vmake
