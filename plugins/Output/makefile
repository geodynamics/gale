#Finds the Absolute path to the Project Root directory
SHELL := /bin/bash
PROJ_ROOT := $(shell until test -r ./Makefile.system ; do cd .. ; done ; echo `pwd`)
include ${PROJ_ROOT}/Makefile.system

# Subdirectories
subdirs := Nusselt $(shell find * -type d -follow -prune ! -name .svn ! -name Nusselt )

include ${PROJ_ROOT}/Makefile.vmake
