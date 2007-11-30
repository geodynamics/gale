
#Finds the Absolute path to the Project Root directory
SHELL := /bin/bash
PROJ_ROOT := $(shell until test -r ./Makefile.system ; do cd .. ; done ; echo `pwd`)
include ${PROJ_ROOT}/Makefile.system

# Subdirectories
tests = lib${PROJECT}_DomainMesh

SRCS = \
	testCartesianGenerator.c \
	testDecomp.c \
	testMesh.c \
	testMeshTopology.c \
	testSpatialTree.c


checks = \
	testCartesianGenerator.0of1.sh \
	testDecomp.0of1.sh \
	testMesh.0of1.sh \
	testCartesianGenerator.0of2.sh \
	testDecomp.0of2.sh \
	testMesh.0of2.sh \
	testCartesianGenerator.0of3.sh \
	testDecomp.0of3.sh \
	testMesh.0of3.sh \
	testCartesianGenerator.0of4.sh \
	testDecomp.0of4.sh \
	testMesh.0of4.sh

EXTERNAL_INCLUDES = -I${INC_DIR}/${PROJECT}
EXTERNAL_LIBS = -L${LIB_DIR} -lStgDomainMesh -lStgDomainShape -lStgDomainGeometry ${RPATH_LFLAGS}

packages = STGERMAIN MPI XML MATH DL 

ifdef HAVE_PYTHON
	packages += PYTHON
endif

ifdef USE_OPTIONAL_MATH
	# For full tensor math it needs some blas lapack functions only if available
	packages := LAPACK_BLAS ${packages}
        ifneq (Darwin,${SYSTEM})
	       EXTERNAL_LIBS +=  -lg2c
	endif
endif

include ${PROJ_ROOT}/Makefile.vmake
