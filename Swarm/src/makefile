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
## $Id: makefile 3462 2006-02-19 06:53:24Z WalterLandry $
##
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# obtain defaults for required variables according to system and project location, and then run the build.
ifndef PROJ_ROOT
	PROJ_ROOT := $(shell until test -r ./Makefile.system ; do cd .. ; done ; echo `pwd`)
endif
include ${PROJ_ROOT}/Makefile.system


lib = StgDomainSwarm
libstatic = yes
libdynamic = yes
includes = StgDomain/Swarm

SRCS = \
	CellLayout.c \
	ElementCellLayout.c \
	FileParticleLayout.c \
	Finalise.c \
	GaussBorderParticleLayout.c \
	GaussParticleLayout.c \
	GlobalParticleLayout.c \
	Init.c \
	IntegrationPoint.c \
	LineParticleLayout.c \
	ManualParticleLayout.c \
	MeshParticleLayout.c \
	OperatorSwarmVariable.c \
	ParticleCommHandler.c \
	ParticleLayout.c \
	ParticleMovementHandler.c \
	ParticleShadowSync.c \
	PerCellParticleLayout.c \
	PlaneParticleLayout.c \
	Random.c \
	RandomParticleLayout.c \
	SingleCellLayout.c \
	SpaceFillerParticleLayout.c \
	StandardParticle.c \
	SwarmClass.c \
	SwarmDump.c \
	SwarmOutput.c \
	SwarmShapeVC.c \
	SwarmVariable.c \
	SwarmVariable_Register.c \
	Swarm_Register.c \
	TriGaussParticleLayout.c \
	TriSingleCellLayout.c \
	UnionParticleLayout.c \
	WithinShapeParticleLayout.c
	
HDRS = \
	CellLayout.h \
	ElementCellLayout.h \
	FileParticleLayout.h \
	Finalise.h \
	GaussBorderParticleLayout.h \
	GaussParticleLayout.h \
	GlobalParticleLayout.h \
	Init.h \
	IntegrationPoint.h \
	LineParticleLayout.h \
	ManualParticleLayout.h \
	MeshParticleLayout.h \
	OperatorSwarmVariable.h \
	ParticleCommHandler.h \
	ParticleLayout.h \
	ParticleMovementHandler.h \
	ParticleShadowSync.h \
	PerCellParticleLayout.h \
	PlaneParticleLayout.h \
	Random.h \
	RandomParticleLayout.h \
	ShadowInfo.h \
	SingleCellLayout.h \
	SpaceFillerParticleLayout.h \
	StandardParticle.h \
	Swarm.h \
	SwarmClass.h \
	SwarmDump.h \
	SwarmOutput.h \
	SwarmShapeVC.h \
	SwarmVariable.h \
	SwarmVariable_Register.h \
	Swarm_Register.h \
	TriGaussParticleLayout.h \
	TriSingleCellLayout.h \
	UnionParticleLayout.h \
	WithinShapeParticleLayout.h \
	shortcuts.h \
	types.h

EXTERNAL_INCLUDES = 
EXTERNAL_LIBS = -L${LIB_DIR} -lStgDomainMesh -lStgDomainShape -lStgDomainGeometry ${RPATH_LFLAGS}

packages = STGERMAIN MPI XML

ifdef HAVE_PYTHON
	packages += PYTHON
endif

include ${PROJ_ROOT}/Makefile.vmake
