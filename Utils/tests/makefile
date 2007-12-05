#Finds the Absolute path to the Project Root directory
SHELL := /bin/bash
PROJ_ROOT := $(shell until test -r ./Makefile.system ; do cd .. ; done ; echo `pwd`)
include ${PROJ_ROOT}/Makefile.system

NAME = $(shell basename `pwd | sed s/tests//g`)

tests = lib${PROJECT}_${NAME}
checks = $(wildcard *.sh)
SRCS = $(wildcard *.c)

EXTERNAL_INCLUDES = -I${INC_DIR}/${PROJECT} 
EXTERNAL_LIBS= -L${LIB_DIR} -l${PROJECT}${NAME} -lStGermainBaseFoundation -lStGermainBaseIO -lStGermainBaseContainer ${RPATH_LFLAGS}

packages = MPI XML MATH DL

include ${PROJ_ROOT}/Makefile.vmake
