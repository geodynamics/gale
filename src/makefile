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

EXTERNAL_LIBS =  -lglucifer ${EXPORT_DYNAMIC_LFLAGS}
EXTERNAL_INCLUDES = -DMODULE_EXT=\"${MODULE_EXT}\" -I${GLUCIFER_INCDIR}

packages = PICELLERATOR STGFEM STGDOMAIN STGERMAIN

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

ifdef HAVE_LIBAVCODEC
	packages += LIBAVCODEC
endif

packages += PETSC MPI XML

include ${PROJ_ROOT}/Makefile.vmake
