	
# obtain defaults for required variables according to system and project location, and then run the build.
override PROJ_ROOT = .
include ${PROJ_ROOT}/Makefile.system

include Makefile.def
subproj = ${def_proj}

include ${PROJ_ROOT}/Makefile.vmake
