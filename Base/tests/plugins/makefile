#Finds the Absolute path to the Project Root directory
SHELL := /bin/bash
PROJ_ROOT := $(shell until test -r ./Makefile.system ; do cd .. ; done ; echo `pwd`)
include ${PROJ_ROOT}/Makefile.system

# External Libraries and Headers
#EXTERNAL_INCLUDES = ${STGERMAIN_INCLUDES} ${MPI_INCLUDES} ${XML_CFLAGS} 
#EXTERNAL_LIBS = ${STGERMAIN_LIBS} ${MPI_LIBS} ${XML_LIBS} 

packages = STGDOMAIN STGERMAIN MPI XML GL

# Note don't put the .c extension. just the base name of the plugin
plugins := $(wildcard *.c)

include ${PROJ_ROOT}/Makefile.vmake

