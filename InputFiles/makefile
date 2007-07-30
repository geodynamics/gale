
SHELL := $(shell which bash)
PROJ_ROOT := $(shell until test -r ./Makefile.system ; do cd .. ; done ; echo `pwd`)

include ${PROJ_ROOT}/Makefile.system

#examples_xml = Arrhenius.xml  Extension.xml  FrankKamenetskii.xml  MobileLid.xml  RayleighTaylorBenchmark.xml ThermoChemBenchmark.xml

modName := $(shell basename `pwd` )

XMLS_TARGET_DIR = ${EXAMPLES_DIR}
XMLS := ${wildcard *.xml}

subdirs := src 

# obtain defaults for required variables according to system and project location, and then run the build.
include ${PROJ_ROOT}/Makefile.vmake
