import os
import config
import config.utils as utils

class ParMETIS(config.Package):

    def __init__(self, ctx, **kw):
        config.Package.__init__(self, ctx, **kw)

    def _setup_dependencies(self):
        config.Package._setup_dependencies(self)
        self.mpi = self.add_dependency(config.packages.MPI,
                                       required=True, combine=True)

    def setup_libraries(self):
        self.add_library_set("parmetis.h", ["parmetis", "metis"])

    source_code = {"c": """#include <stdlib.h>
#include <mpi.h>
#include <parmetis.h>
int main( int argc, char** argv ) {
  return EXIT_SUCCESS;
}
"""}
