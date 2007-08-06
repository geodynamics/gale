
#Finds the Absolute path to the Project Root directory
SHELL := /bin/bash
PROJ_ROOT := $(shell until test -r ./Makefile.system ; do cd .. ; done ; echo `pwd`)
include ${PROJ_ROOT}/Makefile.system

# Subdirectories
subdirs := Base Windowing RenderingEngines OutputFormats InputFormats DrawingObjects WindowInteractions libglucifer src plugins ModelComponents


ifeq (true,$(shell if test -x $(DOXYGEN); then echo true; fi ))
	def_sub += doc
else
	warn := $(shell echo Not entering doc directory as doxygen not found. 1>&2 )
	warn := $(shell echo To generate doxgen docs, please set the DOXYGEN environment variable then re-run configure.sh. 1>&2 )
endif

include ${PROJ_ROOT}/Makefile.vmake
