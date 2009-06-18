import os, glob
import config
import config.utils as utils

class Python(config.Package):

    def __init__(self, ctx, **kw):
        config.Package.__init__(self, ctx, **kw)
        self.version = None

    def setup_libraries(self):
        self.add_library_set(["Python.h"], ["python"])

    def gen_inc_exts(self, inc_dirs):
        ext_map = []
        inc_d = inc_dirs[0]
        for d in glob.glob(os.path.join(inc_d, "python*")):
            if not os.path.isdir(d):
                continue
            d = os.path.basename(d)
            try:
                n = map(int, d[6:].split("."))
                ext_map.append((n, os.path.basename(d)))
            except:
                pass
        if ext_map:
            ext_map.sort(lambda x,y: y[0] < x[0] and -1 or y[0] > x[0] and 1 or 0)
            self.version = ext_map[0][0]
            yield ext_map[0][1]
        for e in config.Package.gen_inc_exts(self, inc_dirs):
            yield e

    def gen_lib_cands(self, dirs):
        if self.version is not None:
            version = ".".join([str(i) for i in self.version])
            for env in self._lib_sets:
                new_libs = []
                for l in env["libs"]:
                    if l == "python":
                        new_libs.append(l + version)
                    else:
                        new_libs.append(l)
                env["libs"] = new_libs
        for c in config.Package.gen_lib_cands(self, dirs):
            yield c

    source_code = {"c": """#include <stdlib.h>
#include <Python.h>
int main(int argc, char** argv) {
  return EXIT_SUCCESS;
}
"""}
