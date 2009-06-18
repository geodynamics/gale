import os
import config
import config.utils as utils


class pcu(config.Package):

    def __init__(self, ctx, **kw):
        config.Package.__init__(self, ctx, **kw)


    def setup_libraries(self):
        self.add_library_set(["pcu/pcu.h"], ["pcu"])


    source_code = {"c": """#include <stdlib.h>
#include <pcu/pcu.h>
int main( int argc, char** argv ) {
  return EXIT_SUCCESS;
}
"""}
