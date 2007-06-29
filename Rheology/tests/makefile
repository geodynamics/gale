
#Finds the Absolute path to the Project Root directory
SHELL := /bin/bash
PROJ_ROOT := $(shell until test -r ./Makefile.system ; do cd .. ; done ; echo `pwd`)
include ${PROJ_ROOT}/Makefile.system

checks = testArrheniusBackground2D.0of1.sh \
	 testByerleeYieldCriterion.0of1.sh \
	testConstitutiveMatrixCartesian2D.0of1.sh \
	testConstitutiveMatrixCartesian3D.0of1.sh \
	testDirector2D_PerMaterial.0of1.sh \
	testDirector2D_RandomDirs.0of1.sh \
	testDirector2D_X.0of1.sh \
	testFrankKamenetskiiBackground2D.0of1.sh \
	testMaterialViscosityBackground2D.0of1.sh \
	testNonNewtonianShearN3.0of1.sh \
	testNonNewtonianShearN4.0of1.sh \
	testSpinUp.0of1.sh \
	testVonMisesYieldCriterion.0of1.sh \
	testLateralViscosityAnalytic2dVisc102.0of1.sh \
	testLateralViscosityAnalytic2dVisc106.0of1.sh \
	#testColumnViscosity.0of1.sh \
	testColumnViscosity.0of2.sh
#$(wildcard *.sh)

subdirs := testConstitutiveMatrix testStiffnessMatrixAssembly testDirector testDirectorRandomness testDirectorPerMaterial testYieldCriterion NonNewtonianShearSolution LateralViscosityAnalytic

include ${PROJ_ROOT}/Makefile.vmake
