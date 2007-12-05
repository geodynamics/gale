# obtain defaults for required variables according to system and project location, and then run the build.
ifndef PROJ_ROOT
	PROJ_ROOT=..
endif
include ${PROJ_ROOT}/Makefile.system

subdirs = src tests

include ${PROJ_ROOT}/Makefile.vmake
