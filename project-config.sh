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
## $Id: project-config.sh 964 2007-10-11 08:03:06Z SteveQuenette $
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

setValue PROJECT 'StgFEM'

# Setup Makefile shortcuts
. ./VMake/Config/makefile-shortcuts.sh

# Do configurations
# Note: StGermain needs to be before MPI and PETSc, due to the
# "compatibility" stuff. -- PatrickSunter
. ./VMake/Config/StGermain-config.sh
. ./VMake/Config/StgDomain-config.sh

. ./VMake/Config/compiler-config.sh
. ./VMake/Config/math-config.sh
. ./VMake/Config/PETSc-config.sh
. ./VMake/Config/mpi-config.sh
. ./VMake/Config/xml-config.sh
. ./VMake/Config/python-config.sh --optional

