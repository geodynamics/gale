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
## $Id: Makefile.rules 2752 2005-03-07 06:45:55Z AlanLo $
##
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#Finds the Absolute path to the Project Root directory
SHELL := /bin/bash
PROJ_ROOT := $(shell until test -r ./Makefile.system ; do cd .. ; done ; echo `pwd`)
include ${PROJ_ROOT}/Makefile.system

NAME = $(shell basename `pwd | sed s/tests//g`)

tests = lib${PROJECT}_${NAME}
SRCS = $(wildcard *.c)

checks = testOutputPPM.0of1.sh
ifdef HAVE_PNG
	checks += testOutputPNG.0of1.sh
endif
ifdef HAVE_JPEG
	checks += testOutputJPEG.0of1.sh
endif
ifdef HAVE_TIFF

ifeq (${SYSTEM_BIG_ENDIAN},0)
	checks += testOutputTIFF-littleendian.0of1.sh
else
	checks += testOutputTIFF-bigendian.0of1.sh
endif

endif
ifdef HAVE_LIBFAME

# The encoded results differ slightly depending on platform
ifneq (${SYSTEM_64},0)
	checks += testEncoderLibfame-64bit.0of1.sh
else
ifeq (${SYSTEM},Darwin)
	checks += testEncoderLibfame-darwin.0of1.sh
else
	checks += testEncoderLibfame.0of1.sh
endif
endif

endif

packages = STGERMAIN MPI XML MATH DL HYPRE

include ${PROJ_ROOT}/Makefile.vmake
