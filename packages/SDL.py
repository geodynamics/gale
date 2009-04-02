import os, platform
import config
import config.utils as utils


class SDL(config.Package):

    def __init__(self, ctx, **kw):
        config.Package.__init__(self, ctx, **kw)
        self.inc_exts = ["SDL"]


    def _setup_dependencies(self):
        config.Package._setup_dependencies(self)
	self.opengl = self.add_dependency(config.packages.OpenGL, required=True, combine=True)
	if platform.system() == "Darwin":
            self.cocoa = self.add_dependency(config.packages.Cocoa, required=False, combine=True)


    def setup_libraries(self):
        self.add_library_set(["SDL.h"], ["SDL"])
        self.add_library_set(["SDL.h"], ["SDLmain", "SDL"])
        self.add_library_set(["SDL.h"], ["SDLmain"], ["SDL"])


    source_code = {"c": """#include <stdlib.h>
#include <SDL.h>
int main( int argc, char** argv ) {
  return EXIT_SUCCESS;
}
"""}
