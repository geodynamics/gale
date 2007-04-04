#Finds the Absolute path to the Project Root directory
SHELL := /bin/bash
PROJ_ROOT := $(shell until test -r ./Makefile.system ; do cd .. ; done ; echo `pwd`)
include ${PROJ_ROOT}/Makefile.system

# Source Code and Header files
SRCS := $(wildcard *.c)
HDRS := $(wildcard *.h)

EXTERNAL_INCLUDES = ${GL_INCLUDES}

# External Libraries and Headers
packages = STGERMAIN MPI XML

# Work out name to call library
lib = $(shell basename `pwd | sed s/src//g`)
libdynamic = yes

# Name of directory to place header files
includes := ${PROJECT}/$(shell basename `pwd | sed s/src//g`)

include ${PROJ_ROOT}/Makefile.vmake
