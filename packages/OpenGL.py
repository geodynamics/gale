import os, platform
import config
import config.utils as utils


class OpenGL(config.Package):

    def __init__(self, ctx):
        config.Package.__init__(self, ctx)
        self.base_dirs = ["/usr/X11R6"]
        self.inc_exts = ["GL"]
	if platform.system() == "Darwin":
            self.base_dirs.append("/System/Library/Frameworks/OpenGL.framework")


    def setup_libraries(self):
        self.add_library_set(["gl.h", "glu.h"], ["GL", "GLU"])
        self.add_library_set(["gl.h", "glu.h"], ["OpenGL"])


    def _test(self, cfg):
        # If we're on Darwin we need this hack.
        if platform.system() == "Darwin":
            cfg.append_unique("lnkprogflags",
                              "-dylib_file /System/Library/Frameworks/" \
                                  "OpenGL.framework/Versions/A/Libraries/libGL.dylib" \
                                  ":/System/Library/Frameworks/OpenGL.framework/Versions" \
                                  "/A/Libraries/libGL.dylib", make_list=True)

        # Now do normal procedure.
        return config.Package._test(self, cfg)


    source_code = {"c": """#include <stdlib.h>
#include <gl.h>
#include <glu.h>
int main( int argc, char** argv ) {
  return EXIT_SUCCESS;
}
"""}
