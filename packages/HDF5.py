import os
import config
import config.utils as utils

class HDF5(config.Package):

    source_code = {"c": """#include <stdlib.h>
#include <hdf5.h>
int main( int argc, char** argv ) {
  return EXIT_SUCCESS;
}
"""}

    def __init__(self, ctx, **kw):
        config.Package.__init__(self, ctx, **kw)
        self.base_patterns = ["*hdf5*", "*HDF5*"]

    def _setup_deps(self):
        config.Package._setup_deps(self)
        self.mpi = self.add_dependency(config.packages.MPI,
                                       required=False, combine=True)

    def setup_libraries(self):
        self.add_library_set("c", ["hdf5.h"], ["hdf5"])
        self.add_auxilliary_libs("c", ["pthread"])
        self.add_auxilliary_libs("c", ["pthread", "z", "sz"])

