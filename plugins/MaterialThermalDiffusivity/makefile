
#Finds the Absolute path to the Project Root directory
SHELL := /bin/bash
PROJ_ROOT := $(shell until test -r ./Makefile.system ; do cd .. ; done ; echo `pwd`)
include ${PROJ_ROOT}/Makefile.system

# Subdirectories
subdirs := $(shell find * -type d -prune ! -name .svn ) 

# Source Code and Header files
SRCS := $(wildcard *.c)
HDRS := $(wildcard *.h)

# What to call the plugin
modName := $(shell basename `pwd | sed s/tests//g`)
mod = ${PROJECT}_${modName}

# Where to put header files
includes = ${PROJECT}/${modName}

# External Libraries and Headers
EXTERNAL_INCLUDES = ${STGFEM_INCLUDES} ${STGERMAIN_INCLUDES} ${PETSC_INCLUDES} ${MPI_INCLUDES} ${XML_CFLAGS} ${PYTHON_INCLUDES} ${PICELLERATOR_INCLUDES}
EXTERNAL_LIBS = \
	${STGFEM_RPATH} ${STGFEM_LIBS} \
	${STGFEM_RPATH} ${STGERMAIN_LIBS} \
	${PETSC_RPATH} ${PETSC_LIBS} \
	${MPI_RPATH} ${MPI_LIBS} \
	${XML_LIBS} ${PYTHON_LIBS}

include ${PROJ_ROOT}/Makefile.vmake
