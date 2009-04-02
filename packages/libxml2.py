import os
import config
import config.utils as utils


class libxml2(config.Package):

    def __init__(self, ctx, **kw):
        config.Package.__init__(self, ctx, **kw)
        self.inc_exts = ["libxml2"]


    def setup_libraries(self):
        inc = os.path.join("libxml", "parser.h")
        self.add_library_set([inc], ["xml2"])


    source_code = {"c": """#include <stdlib.h>
#include <libxml/parser.h>
int main( int argc, char** argv ) {
  return EXIT_SUCCESS;
}
"""}
