import os
import config
import config.utils as utils

class OSMesa(config.Package):

    def __init__(self, ctx, **kw):
        config.Package.__init__(self, ctx, **kw)
        self.base_dirs = ["/usr/X11R6"]
        self.inc_exts = ["GL"]

    def _setup_deps(self):
        config.Package._setup_deps(self)
        self.gl = self.add_dependency(config.packages.OpenGL,
                                       required=True, combine=True)

    def setup_libraries(self):
        self.add_library_set(["osmesa.h"], ["OSMesa"])
        self.add_auxilliary_libs("c", ["dl"])

    source_code = {"c": """#include <stdlib.h>
#include <gl.h>
#include <glu.h>
#include <osmesa.h>
int main( int argc, char** argv ) {
  void* ctx;
  ctx = OSMesaCreateContext(OSMESA_RGBA, NULL);
  OSMesaDestroyContext(ctx);
  return EXIT_SUCCESS;
}
"""}
