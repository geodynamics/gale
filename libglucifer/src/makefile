##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
## Copyright (C), 2003, Victorian Partnership for Advanced Computing (VPAC) Ltd, 110 Victoria Street, Melbourne, 3053, Australia.
##
## Authors:
##	Stevan M. Quenette, Senior Software Engineer, VPAC. (steve@vpac.org)
##	Patrick D. Sunter, Software Engineer, VPAC. (pds@vpac.org)
##	Luke J. Hodkinson, Computational Engineer, VPAC. (lhodkins@vpac.org)
##	Siew-Ching Tan, Software Engineer, VPAC. (siew@vpac.org)
##	Alan H. Lo, Computational Engineer, VPAC. (alan@vpac.org)
##	Raquibul Hassan, Computational Engineer, VPAC. (raq@vpac.org)
##
## This file may be distributed under the terms of the VPAC Public License
## as defined by VPAC of Australia and appearing in the file
## LICENSE.VPL included in the packaging of this file.
##
## This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
## WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
##
## $Id: Makefile.rules 2746 2005-03-06 00:29:00Z SteveQuenette $
##
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#Finds the Absolute path to the Project Root directory
SHELL := /bin/bash
PROJ_ROOT := $(shell until test -r ./Makefile.system ; do cd .. ; done ; echo `pwd`)
include ${PROJ_ROOT}/Makefile.system

lib = ${PROJECT}
libstatic = yes
libdynamic = yes
includes = ${PROJECT}

SRCS = $(wildcard *.c)

HDRS = $(wildcard *.h)

PROJ_OBJLISTS = lucWindowInteractions lucInputFormats lucOutputFormats lucDrawingObjects lucRenderingEngines lucWindowing lucBase
EXTERNAL_LIBS = -L${LIB_DIR} 

EXTERNAL_INCLUDES = -I${INC_DIR}/${PROJECT}  -DLIB_DIR=\"${LIB_DIR}\"

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
ifdef HAVE_TIFF
	packages += TIFF
endif
ifdef HAVE_JPEG
	packages += JPEG
endif
ifdef HAVE_PNG
	packages += PNG
endif
ifdef HAVE_LIBFAME
	packages += LIBFAME
endif
ifdef HAVE_SDL
	packages += SDL
endif
ifdef HAVE_CARBON
	packages += CARBON
endif
ifdef HAVE_LIBAVCODEC
	packages += LIBAVCODEC
endif

ifdef HAVE_VTK
	packages += VTK
endif

packages += PETSC MPI XML MATH 

include ${PROJ_ROOT}/Makefile.vmake
