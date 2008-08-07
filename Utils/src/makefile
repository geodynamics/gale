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


lib = StgDomainUtils
libstatic = yes
libdynamic = yes
includes = StgDomain/Utils

SRCS = \
	AllElementsVC.c \
	AllNodesVC.c \
	CornerVC.c \
	ContactVC.c \
	DomainContext.c \
	DofLayout.c \
	FieldVariable.c \
	FieldVariable_Register.c \
	Finalise.c \
	Init.c \
	InnerWallVC.c \
	LinearRegression.c \
	MeshShapeVC.c \
	MeshBoundaryShape.c \
	NewRemesher.c \
	Operator.c \
	OperatorFieldVariable.c \
	RegularMeshUtils.c \
	RegularRemesher.c \
	RegularRemesherCmpt.c \
	ShapeAdvector.c \
	SobolGenerator.c \
	TimeIntegratee.c \
	TimeIntegrator.c \
	WallVC.c
#	CellRemesher.c \
#	EmbeddedSurface.c \
#	SemiRegDeform.c \
#	StripRemesher.c \
#

HDRS = \
	AllElementsVC.h \
	AllNodesVC.h \
	CornerVC.h \
	ContactVC.h \
	DomainContext.h \
	DofLayout.h \
	FieldVariable.h \
	FieldVariable_Register.h \
	Finalise.h \
	Init.h \
	InnerWallVC.h \
	LinearRegression.h \
	MeshShapeVC.h \
	MeshBoundaryShape.h \
	NewRemesher.def \
	NewRemesher.h \
	Operator.h \
	OperatorFieldVariable.h \
	RegularMeshUtils.h \
	RegularRemesher.def \
	RegularRemesher.h \
	RegularRemesherCmpt.h \
	ShapeAdvector.h \
	SobolGenerator.h \
	TimeIntegratee.h \
	TimeIntegrator.h \
	Utils.h \
	WallVC.h \
	types.h
#	CellRemesher.h \
#	EmbeddedSurface.h \
#	SemiRegDeform.h \
#	StripRemesher.h \
#


EXTERNAL_INCLUDES = 
EXTERNAL_LIBS = -L${LIB_DIR} -lStgDomainMesh -lStgDomainShape -lStgDomainGeometry ${RPATH_LFLAGS}

packages = STGERMAIN MPI XML

ifdef HAVE_PYTHON
	packages += PYTHON
endif

include ${PROJ_ROOT}/Makefile.vmake
