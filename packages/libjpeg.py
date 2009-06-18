import os
import config
import config.utils as utils

class libjpeg(config.Package):

    source_code = {"c": """#include <stdlib.h>
#include <stdio.h>
#include <jpeglib.h>
int main( int argc, char** argv ) {
  return EXIT_SUCCESS;
}
"""}

    def __init__(self, ctx, **kw):
        config.Package.__init__(self, ctx, **kw)

    def setup_libraries(self):
        self.add_library_set(["jpeglib.h"], ["jpeg"])
