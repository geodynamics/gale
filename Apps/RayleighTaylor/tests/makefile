##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
## Copyright (C), 2003-2006, Victorian Partnership for Advanced Computing (VPAC) Ltd, 110 Victoria Street,
##	Melbourne, 3053, Australia.
## Copyright (c) 2005-2006, Monash Cluster Computing, Building 28, Monash University Clayton Campus,
##	Victoria, 3800, Australia
##
## Primary Contributing Organisations:
##	Victorian Partnership for Advanced Computing Ltd, Computational Software Development - http://csd.vpac.org
##	Australian Computational Earth Systems Simulator - http://www.access.edu.au
##	Monash Cluster Computing - http://www.mcc.monash.edu.au
##
## Contributors:
##	Robert Turnbull, Research Assistant, Monash University. (robert.turnbull@sci.monash.edu.au)
##	Patrick D. Sunter, Software Engineer, VPAC. (patrick@vpac.org)
##	Alan H. Lo, Computational Engineer, VPAC. (alan@vpac.org)
##	Stevan M. Quenette, Senior Software Engineer, VPAC. (steve@vpac.org)
##	David May, PhD Student, Monash University (david.may@sci.monash.edu.au)
##	Vincent Lemiale, Postdoctoral Fellow, Monash University. (vincent.lemiale@sci.monash.edu.au)
##	Julian Giordani, Research Assistant, Monash University. (julian.giordani@sci.monash.edu.au)
##	Louis Moresi, Associate Professor, Monash University. (louis.moresi@sci.monash.edu.au)
##	Luke J. Hodkinson, Computational Engineer, VPAC. (lhodkins@vpac.org)
##	Raquibul Hassan, Computational Engineer, VPAC. (raq@vpac.org)
##	David Stegman, Postdoctoral Fellow, Monash University. (david.stegman@sci.monash.edu.au)
##	Wendy Sharples, PhD Student, Monash University (wendy.sharples@sci.monash.edu.au)
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
## $Id: makefile 409 2007-01-17 13:01:26Z SteveQuenette $
##
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# obtain defaults for required variables according to system and project location, and then run the build.
SHELL := $(shell which bash)
PROJ_ROOT := $(shell until test -r ./Makefile.system ; do cd .. ; done ; echo `pwd`)

include ${PROJ_ROOT}/Makefile.system

subdirs = AnalyticColumn AnalyticBlock

modName = $(shell basename `pwd | sed s/tests//g`)

checks = \
	testAnalyticBlock.0of1.sh \
	testAnalyticBlock.0of2.sh \
	testAnalyticColumn.0of1.sh \
	testAnalyticColumn.0of2.sh \
	testAnalyticColumn3D.0of1.sh \
	testRayTay3D.0of1.sh \
	testRayTaySuperimposedBall2D.0of1.sh \
	testRayTaySuperimposedBall2D-invert.0of1.sh

ifeq (${MACHINE},i686)
	# i686 archs have weaker floating point calculations
	checks += testRayTay2D-i686.0of1.sh 
else
	checks += testRayTay2D.0of1.sh 
endif

EXTERNAL_INCLUDES = -I${INC_DIR}/${PROJECT} 

XMLS_TARGET_DIR = ${EXAMPLES_DIR}/${modName}/tests
XMLS = $(wildcard *.xml)

packages = STGERMAIN PETSC MPI XML MATH DL HYPRE

include ${PROJ_ROOT}/Makefile.vmake
