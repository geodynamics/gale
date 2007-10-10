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

# obtain defaults for required variables according to system and project location, and then run the build.
ifndef PROJ_ROOT
	PROJ_ROOT := $(shell until test -r ./Makefile.system ; do cd .. ; done ; echo `pwd`)
endif
include ${PROJ_ROOT}/Makefile.system


lib = StgDomainShape
libstatic = yes
libdynamic = yes
includes = StgDomain/Shape

SRCS = \
	BelowCosinePlane.c \
	Cylinder.c \
	Intersection.c \
	Sphere.c \
	BelowPlane.c \
	Everywhere.c \
	PolygonShape.c \
	Superellipsoid.c \
	Box.c \
	Finalise.c \
	PythonShape.c \
	Union.c \
	ConvexHull.c \
	Init.c \
	ShapeClass.c


HDRS = \
	BelowCosinePlane.h \
	Cylinder.h \
	Intersection.h \
	ShapeClass.h \
	types.h \
	BelowPlane.h \
	Everywhere.h \
	PolygonShape.h \
	Sphere.h \
	Box.h \
	Finalise.h \
	PythonShape.h \
	Superellipsoid.h \
	ConvexHull.h \
	Init.h \
	Shape.h \
	Union.h

EXTERNAL_INCLUDES = 
EXTERNAL_LIBS = -L${LIB_DIR} -lStgDomainGeometry ${RPATH_LFLAGS}

packages = STGERMAIN MPI XML

ifdef HAVE_PYTHON
	packages += PYTHON
endif

include ${PROJ_ROOT}/Makefile.vmake
