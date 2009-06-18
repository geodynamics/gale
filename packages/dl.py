import os
import config
import config.utils as utils

class dl(config.Package):

    source_code = {"c": """#include <stdlib.h>
#include <dlfcn.h>
int main( int argc, char** argv ) {
  return EXIT_SUCCESS;
}
"""}

    def __init__(self, ctx, **kw):
        config.Package.__init__(self, ctx, **kw)

    def setup_libraries(self):
        self.add_library_set("c", ["dlfcn.h"], ["dl"])
