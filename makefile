	
# obtain defaults for required variables according to system and project location, and then run the build.
override PROJ_ROOT = .
include ${PROJ_ROOT}/Makefile.system

subproj = StGermain StgDomain StgFEM PICellerator gLucifer Underworld 

include ${PROJ_ROOT}/Makefile.vmake
