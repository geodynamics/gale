import os
import config

class PICellerator(config.Package):

    def __init__(self, ctx, **kw):
        config.Package.__init__(self, ctx, **kw)
        self.inc_exts = ["PICellerator"]

    def setup_dependencies(self):
        config.Package.setup_dependencies(self)
        self.fem = self.add_dependency(config.packages.StgFEM, required=True, combine=True)

    def setup_libraries(self):
        self.add_library_set([os.path.join("PICellerator", "PICellerator.h")], ["PICellerator"])

    def process_options(self):
        config.Package.process_options(self)
        self.stg = self.fem.dom.stg
        base_dir = self.stg.get_option("stg_dir", None)
        if base_dir is not None:
            self.forced_base_dirs = [base_dir]

    source_code = {"c": """#include <stdlib.h>
#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>
#include <StgFEM/StgFEM.h>
#include <PICellerator/PICellerator.h>
int main( int argc, char** argv ) {
  return EXIT_SUCCESS;
}
"""}
