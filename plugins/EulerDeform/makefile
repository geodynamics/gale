SHELL := $(shell which bash)
PROJ_ROOT := $(shell until test -r ./Makefile.system ; do cd .. ; done ; echo `pwd`)

include ${PROJ_ROOT}/Makefile.system

modName := $(shell basename `pwd` )
modName := $(shell echo ${modName} | cut -c 1 | tr '[:lower:]' '[:upper:]' )$(shell echo ${modName} | cut -c 2- ) 

mod = ${PROJECT}_${modName}
includes = ${PROJECT}/${modName}

SRCS = $(wildcard *.c)
HDRS = $(wildcard *.h)

packages = PICELLERATOR STGFEM STGDOMAIN STGERMAIN PETSC MPI XML

include ${PROJ_ROOT}/Makefile.vmake
