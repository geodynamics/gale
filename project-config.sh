#!/bin/sh
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
## Victorian Partnership for Advanced Computing (VPAC) Ltd, Australia
## (C) 2003 All Rights Reserved
##
## California Institute of Technology (Caltech), USA
## (C) 2003 All Rights Reserved
##
## Authors:
## 	Stevan M. Quenette, Senior Software Engineer, VPAC.
##	Stevan M. Quenette, Visitor in Geophysics, Caltech.
##
## <copyright-release-tag>
##
## Role:
##	Obtain the project configuration.
##
## $Id: project-config.sh 376 2006-10-18 06:58:41Z SteveQuenette $
##
##
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Be Bourne compatible
if test -n "${ZSH_VERSION+set}" && (emulate sh) >/dev/null 2>&1; then
	emulate sh
	NULLCMD=:
elif test -n "${BASH_VERSION+set}" && (set -o posix) >/dev/null 2>&1; then
	set -o posix
fi

. './build-functions.sh'

setValue PROJECT 'PICellerator'

# Setup Makefile shortcuts
. ./VMake/Config/makefile-shortcuts.sh

# Do configurations
. ./VMake/Config/StGermain-config.sh
. ./VMake/Config/StgFEM-config.sh

. ./VMake/Config/compiler-config.sh
. ./VMake/Config/math-config.sh
. ./VMake/Config/PETSc-config.sh
. ./VMake/Config/mpi-config.sh
. ./VMake/Config/xml-config.sh
. ./VMake/Config/python-config.sh --optional
