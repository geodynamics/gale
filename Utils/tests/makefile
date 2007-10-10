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
##  This library is free software; you can redistribute it and/or
##  modify it under the terms of the GNU Lesser General Public
##  License as published by the Free Software Foundation; either
##  version 2.1 of the License, or (at your option) any later version.
##
##  This library is distributed in the hope that it will be useful,
##  but WITHOUT ANY WARRANTY; without even the implied warranty of
##  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
##  Lesser General Public License for more details.
##
##  You should have received a copy of the GNU Lesser General Public
##  License along with this library; if not, write to the Free Software
##  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
##
## $Id: makefile 3982 2007-01-23 03:52:15Z KathleenHumble $
##
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#Finds the Absolute path to the Project Root directory
SHELL := /bin/bash
	PROJ_ROOT := $(shell until test -r ./Makefile.system ; do cd .. ; done ; echo `pwd`)
include ${PROJ_ROOT}/Makefile.system

tests = lib${PROJECT}_Utils

SRCS = \
	testRegularMeshUtils.c \
	testFieldVariable_Register.c \
	testOperatorFieldVariable.c \
	testAllNodesVC.c \
	testWallVC.c \
	testCompositeVC.c \
	testCompositeVC-dictionary.c \
	testDofLayout.c \
	testDofLayout-remap.c \
	testDofLayout-saveAndLoad.c \
	testSobolGenerator.c \
	testTimeIntegration.c \
	testCornerVC.c \
	testRegularRemesher.c
#	testSemiRegDeform.c \

checks = \
	testFieldVariable_Register.0of1.sh \
	testOperatorFieldVariable.0of1.sh \
	testAllNodesVC.0of1.sh \
	testAllNodesVC.0of2.sh \
	testAllNodesVC.1of2.sh \
	testAllNodesVC.0of3.sh \
	testAllNodesVC.1of3.sh \
	testAllNodesVC.2of3.sh \
	testWallVC.0of1.sh \
	testWallVC.0of2.sh \
	testWallVC.1of2.sh \
	testWallVC.0of3.sh \
	testWallVC.1of3.sh \
	testWallVC.2of3.sh \
	testCompositeVC.0of1.sh \
	testCompositeVC.0of2.sh \
	testCompositeVC.1of2.sh \
	testCompositeVC.0of3.sh \
	testCompositeVC.1of3.sh \
	testCompositeVC.2of3.sh \
	testCompositeVC-dictionary.0of1.sh \
	testDofLayout.0of1.sh \
	testDofLayout-remap.0of1.sh \
	testDofLayout-saveAndLoad.0of1.sh \
	testTimeIntegrationEuler.0of1.sh \
	testTimeIntegrationRK2.0of1.sh \
	testTimeIntegrationRK4.0of1.sh \
	testSobolGenerator.0of1.sh \
	testRegularMeshUtils.0of1.sh \
	testRegularMeshUtils.0of2.sh \
	testRegularMeshUtils.1of2.sh \
	testCornerVC.0of1.sh \
	testCornerVC.0of2.sh \
	testCornerVC.1of2.sh \
	testCornerVC.0of3.sh \
	testCornerVC.1of3.sh \
	testCornerVC.2of3.sh

PROJ_LIBS = 
EXTERNAL_INCLUDES = -I${INC_DIR}/${PROJECT}
EXTERNAL_LIBS = -L${LIB_DIR} -lStgDomainUtils -lStgDomainMesh -lStgDomainShape -lStgDomainGeometry ${RPATH_LFLAGS}

packages = STGERMAIN MPI XML MATH DL 

ifdef HAVE_PYTHON
	packages += PYTHON
endif

ifdef USE_OPTIONAL_MATH
	# For full tensor math it needs some blas lapack functions only if available
	packages := LAPACK_BLAS ${packages}
        ifneq (Darwin,${SYSTEM})
	       EXTERNAL_LIBS +=  -lg2c
	endif
endif

include ${PROJ_ROOT}/Makefile.vmake
