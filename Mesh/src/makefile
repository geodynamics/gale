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
	PROJ_ROOT=../..
endif
include ${PROJ_ROOT}/Makefile.system


lib = StgDomainMesh
libstatic = yes
libdynamic = yes
includes = StgDomain/Mesh

SRCS = 				\
	Grid.c				\
	Decomp.c			\
	Sync.c				\
	MeshTopology.c			\
	IGraph.c			\
	Mesh_ElementType.c		\
	Mesh_HexType.c			\
	Mesh_CentroidType.c		\
	Mesh_Algorithms.c		\
	Mesh_HexAlgorithms.c		\
	Mesh_RegularAlgorithms.c	\
	Mesh_CentroidAlgorithms.c	\
	MeshClass.c			\
	MeshGenerator.c			\
	CartesianGenerator.c		\
	MeshAdaptor.c			\
	SurfaceAdaptor.c		\
	CompressionAdaptor.c		\
	MeshVariable.c			\
	Init.c				\
	Finalise.c

HDRS = 				\
	types.h				\
	shortcuts.h			\
	Grid.h				\
	Decomp.def			\
	Decomp.h			\
	Sync.def			\
	Sync.h				\
	MeshTopology.def		\
	MeshTopology.h			\
	IGraph.def			\
	IGraph.h			\
	Mesh.h				\
	Mesh_ElementType.h		\
	Mesh_HexType.h			\
	Mesh_CentroidType.h		\
	Mesh_Algorithms.h		\
	Mesh_HexAlgorithms.h		\
	Mesh_RegularAlgorithms.h	\
	Mesh_CentroidAlgorithms.h	\
	MeshClass.h			\
	MeshGenerator.h			\
	CartesianGenerator.h		\
	MeshAdaptor.h			\
	SurfaceAdaptor.h		\
	CompressionAdaptor.h		\
	MeshVariable.h			\
	Init.h				\
	Finalise.h

EXTERNAL_INCLUDES = 
EXTERNAL_LIBS = -L${LIB_DIR} -lStgDomainShape -lStgDomainGeometry ${RPATH_LFLAGS}

packages = STGERMAIN MPI XML

ifdef HAVE_PYTHON
	packages += PYTHON
endif

include ${PROJ_ROOT}/Makefile.vmake
