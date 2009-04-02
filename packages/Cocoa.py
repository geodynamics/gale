import os, platform
import config
import config.utils as utils


class Cocoa(config.Package):

    def __init__(self, ctx):
        config.Package.__init__(self, ctx)
        self.base_dirs.append("/System/Library/Frameworks/Cocoa.framework")


    def setup_libraries(self):
        self.add_library_set(["Cocoa.h"], ["Cocoa"])


    source_code = {"c": """#include <stdlib.h>
int main( int argc, char** argv ) {
  return EXIT_SUCCESS;
}
"""}
