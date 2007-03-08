SHELL := $(shell which bash)
PROJ_ROOT := $(shell until test -r ./Makefile.system ; do cd .. ; done ; echo `pwd`)

include ${PROJ_ROOT}/Makefile.system

#subdirs = tests

modName := $(shell basename `pwd` )
mod = ${PROJECT}_${modName}
includes = ${PROJECT}/${modName}

SRCS = $(wildcard *.c)
HDRS = $(wildcard *.h)

packages = PICELLERATOR STGFEM STGERMAIN PETSC MPI XML PYTHON

# obtain defaults for required variables according to system and project location, and then run the build.
include ${PROJ_ROOT}/Makefile.vmake
