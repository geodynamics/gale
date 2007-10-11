#Finds the Absolute path to the Project Root directory
SHELL := /bin/bash
PROJ_ROOT := $(shell until test -r ./Makefile.system ; do cd .. ; done ; echo `pwd`)
include ${PROJ_ROOT}/Makefile.system

# Source Code and Header files
SRCS := $(wildcard *.c)
HDRS := $(wildcard *.h)

packages = STGDOMAIN STGERMAIN

ifdef HAVE_OSMESA
	packages += OSMESA
endif

ifdef HAVE_GL
	packages += GL
endif

packages += MPI XML 

# Work out name to call library
lib = luc$(shell basename `pwd | sed s/src//g`)
libdynamic = yes
libstatic = yes

# Name of directory to place header files
includes := ${PROJECT}/$(shell basename `pwd | sed s/src//g`)

include ${PROJ_ROOT}/Makefile.vmake
