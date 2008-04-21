
#Finds the Absolute path to the Project Root directory
SHELL := /bin/bash
PROJ_ROOT := $(shell until test -r ./Makefile.system ; do cd .. ; done ; echo `pwd`)
include ${PROJ_ROOT}/Makefile.system

# External Libraries and Headers
EXTERNAL_INCLUDES = -I${INC_DIR}/${PROJECT} ${GL_INCLUDES}
EXTERNAL_LIBS = ${RPATH_LFLAGS} -L${LIB_DIR} -llucBase -llucRenderingEngines -lDummyOpenGL

subdirs = DummyOpenGL

packages = STGDOMAIN STGERMAIN GL PETSC MPI XML 

ifdef HAVE_GL2PS
	packages += GL2PS
endif

tests = lib${PROJECT}_RenderingEngines
SRCS = $(wildcard *.c)
checks = $(wildcard *.sh)

include ${PROJ_ROOT}/Makefile.vmake
