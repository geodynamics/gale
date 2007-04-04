#Finds the Absolute path to the Project Root directory
SHELL := /bin/bash
PROJ_ROOT := $(shell until test -r ./Makefile.system ; do cd .. ; done ; echo `pwd`)
include ${PROJ_ROOT}/Makefile.system

# Source Code and Header files
SRCS := $(wildcard *.c) $(wildcard *.cpp)
HDRS := $(wildcard *.h)

# External Libraries and Headers
EXTERNAL_INCLUDES = -I${INC_DIR}
EXTERNAL_LIBS = -L${LIB_DIR} -lglucifer ${GLUCIFER_RPATH} 

packages = PICELLERATOR STGFEM STGERMAIN

ifdef HAVE_OSMESA
	packages += OSMESA
endif
ifdef HAVE_GL
	packages += GL
endif
ifdef HAVE_X11
	packages += X11
endif
ifdef HAVE_VTK
	packages += VTK
endif

packages += PETSC MPI XML

# Work out name to call library
mod = $(shell basename `pwd | sed s/src//g`)

# Name of directory to place header files
includes := ${PROJECT}/${mod}

include ${PROJ_ROOT}/Makefile.vmake
