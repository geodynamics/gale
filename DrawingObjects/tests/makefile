
#Finds the Absolute path to the Project Root directory
SHELL := /bin/bash
PROJ_ROOT := $(shell until test -r ./Makefile.system ; do cd .. ; done ; echo `pwd`)
include ${PROJ_ROOT}/Makefile.system

# External Libraries and Headers
EXTERNAL_INCLUDES = -I${INC_DIR}/${PROJECT} ${GL_INCLUDES}
EXTERNAL_LIBS = ${RPATH_LFLAGS} -L${LIB_DIR} -llucBase -llucWindowing -llucRenderingEngines  -llucOutputFormats -llucInputFormats -llucDrawingObjects -lDummyOpenGL

subdirs = DummyFieldVariable DummySwarmVariable 

packages = PICELLERATOR STGFEM STGERMAIN 

ifdef HAVE_OSMESA
	packages += OSMESA
endif
ifdef HAVE_GL
	packages += GL
endif

packages += PETSC MPI XML 

tests = lib${PROJECT}_DrawingObjects
SRCS = $(wildcard *.c)
checks = $(wildcard test*.sh)

include ${PROJ_ROOT}/Makefile.vmake
