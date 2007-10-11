##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
## Copyright (C), 2003-2006, Victorian Partnership for Advanced Computing (VPAC) Ltd, 110 Victoria Street,
##	Melbourne, 3053, Australia.
##
## Primary Contributing Organisations:
##	Victorian Partnership for Advanced Computing Ltd, Computational Software Development - http://csd.vpac.org
##	Australian Computational Earth Systems Simulator - http://www.access.edu.au
##	Monash Cluster Computing - http://www.mcc.monash.edu.au
##	Computational Infrastructure for Geodynamics - http://www.geodynamics.org
##
## Contributors:
##	Patrick D. Sunter, Software Engineer, VPAC. (pds@vpac.org)
##	Robert Turnbull, Research Assistant, Monash University. (robert.turnbull@sci.monash.edu.au)
##	Stevan M. Quenette, Senior Software Engineer, VPAC. (steve@vpac.org)
##	David May, PhD Student, Monash University (david.may@sci.monash.edu.au)
##	Louis Moresi, Associate Professor, Monash University. (louis.moresi@sci.monash.edu.au)
##	Luke J. Hodkinson, Computational Engineer, VPAC. (lhodkins@vpac.org)
##	Alan H. Lo, Computational Engineer, VPAC. (alan@vpac.org)
##	Raquibul Hassan, Computational Engineer, VPAC. (raq@vpac.org)
##	Julian Giordani, Research Assistant, Monash University. (julian.giordani@sci.monash.edu.au)
##	Vincent Lemiale, Postdoctoral Fellow, Monash University. (vincent.lemiale@sci.monash.edu.au)
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
## $Id: makefile 964 2007-10-11 08:03:06Z SteveQuenette $
##
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# obtain defaults for required variables according to system and project location, and then run the build.
SHELL := $(shell which bash)
PROJ_ROOT := $(shell until test -r ./Makefile.system ; do cd .. ; done ; echo `pwd`)

include ${PROJ_ROOT}/Makefile.system

subdirs = LidDrivenIsoviscousAnalytic SimpleShearAnalytic
tests = ${PROJECT}-StokesMomentumUzawaTest

checks = \
	testStokesMomentum-UzawaSolver-2D.0of1.sh \
	testStokesMomentum-UzawaSolver-2D.0of2.sh \
	testStokesMomentum-UzawaSolver-3D.0of1.sh \
	testStokesMomentum-UzawaSolver-3D.0of2.sh \
	testStokesMomentum-UzawaSolver-3D.0of3.sh \
	testStokesMomentum-UzawaSolver-SimpleShear2D.0of1.sh \
	testStokesMomentum-UzawaSolver-SimpleShear2D.0of2.sh \
	testStokesMomentum-UzawaSolver-SimpleShear2D.0of3.sh \
	testStokesMomentum-UzawaSolver-SimpleShear3D.0of1.sh \
	testStokesMomentum-UzawaSolver-SimpleShear3D.0of2.sh

ifeq (${PETSC_VERSION_RELEASE}, 0)
	checks += \
		testStokesMomentum-UzawaSolver-2D-MG.0of1.sh \
		testStokesMomentum-UzawaSolver-2D-MG.0of2.sh
endif


PROJ_LIBS =
EXTERNAL_INCLUDES = -I${INC_DIR}/${PROJECT} 

packages = STGDOMAIN STGERMAIN PETSC MPI XML MATH DL HYPRE

include ${PROJ_ROOT}/Makefile.vmake
