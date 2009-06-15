import os
import config
import config.utils as utils

class X11(config.Package):

    def __init__(self, ctx, **kw):
        config.Package.__init__(self, ctx, **kw)
        self.base_dirs = ["/usr/X11R6"]
        self.inc_exts = ["X11"]

    def _setup_dependencies(self):
        config.Package._setup_dependencies(self)
        self.gl = self.add_dependency(config.packages.OpenGL, required=True, combine=True)

    def setup_libraries(self):
        self.add_library_set(["Xlib.h"], ["X11", "Xmu", "xcb", "xcb-xlib", "Xau", "Xdmcp"])
        self.add_library_set(["Xlib.h"], ["X11", "Xmu", "Xau", "Xdmcp"])


    source_code = {"c": """#include <stdlib.h>
#include <Xlib.h>
int main( int argc, char** argv ) {
  void* display;
  display = XOpenDisplay(NULL);
  return EXIT_SUCCESS;
}
"""}
