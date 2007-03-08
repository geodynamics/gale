
SHELL := $(shell which bash)
PROJ_ROOT := $(shell until test -r ./Makefile.system ; do cd .. ; done ; echo `pwd`)

include ${PROJ_ROOT}/Makefile.system

XMLS_TARGET_DIR = ${EXPORT_INCLUDE_DIR}
XMLS := ${wildcard *.xml}

subdirs := VariableConditions BaseApps Viewports 

# obtain defaults for required variables according to system and project location, and then run the build.
include ${PROJ_ROOT}/Makefile.vmake
