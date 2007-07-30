
override PROJ_ROOT=.
include ${PROJ_ROOT}/Makefile.system

subdirs = Rheology Utils libUnderworld src plugins InputFiles SysTest

ifeq (true,$(shell if test -x $(DOXYGEN); then echo true; fi ))
	subdirs += doc
else
	warn := $(shell echo Not entering doc directory as doxygen not found. 1>&2 )
	warn := $(shell echo To generate doxgen docs, please set the DOXYGEN environment variable then re-run configure.sh. 1>&2 )
endif

include ${PROJ_ROOT}/Makefile.vmake


