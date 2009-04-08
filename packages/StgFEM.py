import os
import config

class StgFEM(config.Package):

    def __init__(self, ctx, **kw):
        config.Package.__init__(self, ctx, **kw)
        self.inc_exts = ["StgFEM"]

    def setup_dependencies(self):
        config.Package.setup_dependencies(self)
        self.add_dependency(config.packages.PETSc, required=True, combine=True)
        self.dom = self.add_dependency(config.packages.StgDomain, required=True, combine=True)

    def setup_libraries(self):
        self.add_library_set([os.path.join("StgFEM", "StgFEM.h")], ["StgFEM"])

    def process_options(self):
        config.Package.process_options(self)
        self.stg = self.dom.stg
        base_dir = self.stg.get_option("stg_dir", None)
        if base_dir is not None:
            self.forced_base_dirs = [base_dir]

    source_code = {"c": """#include <stdlib.h>
#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>
#include <StgFEM/StgFEM.h>
int main( int argc, char** argv ) {
  return EXIT_SUCCESS;
}
"""}
