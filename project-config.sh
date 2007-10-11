#!/bin/sh

# Be Bourne compatible
if test -n "${ZSH_VERSION+set}" && (emulate sh) >/dev/null 2>&1; then
	emulate sh
	NULLCMD=:
elif test -n "${BASH_VERSION+set}" && (set -o posix) >/dev/null 2>&1; then
	set -o posix
fi

. './build-functions.sh'

setValue  PROJECT 'Underworld'

# Setup Makefile shortcuts
. ./VMake/Config/makefile-shortcuts.sh

# Do configurations
. ./VMake/Config/compiler-config.sh
. ./VMake/Config/PETSc-config.sh
. ./VMake/Config/mpi-config.sh
. ./VMake/Config/xml-config.sh

. ./VMake/Config/StGermain-config.sh
. ./VMake/Config/StgDomain-config.sh
. ./VMake/Config/StgFEM-config.sh
. ./VMake/Config/PICellerator-config.sh
. ./VMake/Config/gLucifer-config.sh
. ./VMake/Config/sdl-config.sh
. ./VMake/Config/SLEPC-config.sh --optional
