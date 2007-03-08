
#Finds the Absolute path to the Project Root directory
SHELL := /bin/bash
PROJ_ROOT := $(shell until test -r ./Makefile.system ; do cd .. ; done ; echo `pwd`)
include ${PROJ_ROOT}/Makefile.system

checks = $(wildcard *.sh)

subdirs := testConstitutiveMatrix testStiffnessMatrixAssembly LateralViscosityAnalytic  testDirector testDirectorRandomness testDirectorPerMaterial testYieldCriterion NonNewtonianShearSolution

include ${PROJ_ROOT}/Makefile.vmake
