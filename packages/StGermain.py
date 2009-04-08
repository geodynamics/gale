import os
import config

class StGermain(config.Package):

    def __init__(self, ctx, **kw):
        config.Package.__init__(self, ctx, **kw)
        self.inc_exts = ["StGermain"]

    def setup_dependencies(self):
        config.Package.setup_dependencies(self)
        self.add_dependency(config.packages.MPI, required=True, combine=True)
        self.add_dependency(config.packages.libxml2, required=True, combine=True)

    def setup_libraries(self):
        self.add_library_set([os.path.join("StGermain", "StGermain.h")], ["StGermain"])

    def setup_options(self):
        config.Package.setup_options(self)
        self.options.add_option(
            config.utils.options.Option("stg_dir", "--stg-dir",
                                        pref_sep="=", help=self.help["stg_dir"]))

    def process_options(self):
        config.Package.process_options(self)
        base_dir = self.get_option("stg_dir", None)
        if base_dir is not None:
            self.forced_base_dirs = [base_dir]

    def make_help_dict(self):
        config.Package.make_help_dict(self)
        self.help["stg_dir"] = "Specify the base location for all StG* related projects."

    source_code = {"c": """#include <stdlib.h>
#include <StGermain/StGermain.h>
int main( int argc, char** argv ) {
  return EXIT_SUCCESS;
}
"""}
