##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
## Victorian Partnership for Advanced Computing (VPAC) Ltd, Australia
## (C) 2003 All Rights Reserved
##
## <copyright-release-tag>
##
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

SHELL := /bin/bash
PROJ_ROOT := $(shell until test -r ./Makefile.system ; do cd .. ; done ; echo `pwd`)
include ${PROJ_ROOT}/Makefile.system

bin = ${PROJECT}
SRCS = main.c

EXTERNAL_LIBS = ${EXPORT_DYNAMIC_LFLAGS}
EXTERNAL_INCLUDES = -DMODULE_EXT=\"${MODULE_EXT}\"

packages = STGERMAIN PETSC MPI XML PYTHON MATH DL

include ${PROJ_ROOT}/Makefile.vmake
