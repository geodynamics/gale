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

subdirs = SingleAttractor VaryingCornerAttractors Bouncer Spinner

tests = lib${PROJECT}_Swarm

SRCS = \
	testElementCellLayout0.c \
	testGaussLayout.c \
	testGaussLayoutSingleCell.c \
	testManualParticleLayout.c \
	testParticleCoords.c \
	testSingleCellLayout0.c \
	testSpaceFillerParticleLayout.c \
	testSwarm.c \
	testSwarmDumpAndLoad.c \
	testSwarmOutput.c \
	testSwarmParticleAdvection.c \
	testSwarmParticleShadowSync.c \
	testWithinShapeParticleLayout.c \

checks = \
	testElementCellLayout0.0of1.sh \
	testElementCellLayout0.0of2.sh \
	testElementCellLayout0.1of2.sh \
	testGaussLayout.0of1.sh \
	testGaussLayoutSingleCell.0of1.sh \
	testLineParticleLayout.0of1.sh \
	testManualParticleLayout.0of1.sh \
	testManualParticleLayout.0of2.sh \
	testManualParticleLayout.1of2.sh \
	testSingleCellLayout0.0of1.sh \
	testSpaceFillerParticleLayout.0of1.sh \
	testSpaceFillerParticleLayout.0of2.sh \
	testSpaceFillerParticleLayout.1of2.sh \
	testSwarm.0of1.sh \
	testSwarm.0of2.sh \
	testSwarm.1of2.sh \
	testSwarmDumpAndLoad.0of1.sh \
	testSwarmDumpAndLoad.0of2.sh \
	testSwarmDumpAndLoad.1of2.sh \
	testSwarmOutput.0of1.sh \
	testSwarmOutput.0of2.sh \
	testSwarmParticleAdvectionSingleAttractor.0of1.sh \
	testSwarmParticleAdvectionSingleAttractor.0of2.sh \
	testSwarmParticleAdvectionSingleAttractor.0of3.sh \
	testSwarmParticleAdvectionSingleAttractor.1of2.sh \
	testSwarmParticleAdvectionSingleAttractor.1of3.sh \
	testSwarmParticleAdvectionSingleAttractor.2of3.sh \
	testSwarmParticleShadowSync.0of2.sh \
	testSwarmParticleShadowSync.0of3.sh \
	testWithinShapeParticleLayout.0of1.sh \
	testWithinShapeParticleLayout.0of2.sh \
	testWithinShapeParticleLayout.1of2.sh \

PROJ_LIBS = 
EXTERNAL_INCLUDES = -I${INC_DIR}/${PROJECT} -DLIB_DIR=\"${LIB_DIR}\"
EXTERNAL_LIBS = -L${LIB_DIR} -lStgDomainSwarm -lStgDomainUtils -lStgDomainMesh -lStgDomainShape -lStgDomainGeometry ${RPATH_LFLAGS}

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
