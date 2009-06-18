import os
import config
import config.utils as utils
from config.utils.options import Option

class MPICompiler(config.Tool):

    def __init__(self, ctx):
        config.Tool.__init__(self, ctx)

    @staticmethod
    def same_candidates(c0, c1):
        if "command" not in c1:
            return False
        return c0["command"] == c1["command"]

    def _test(self, cfg):
        if not self.test_options(cfg):
            return False

        # Get hold of the libraries.
        res, out, err = self.execute(cfg, args=[("mpi_show", True)])
        opts = cfg["options"].parse_option_string(out)
        opts = cfg["options"].gather(opts)

        # Remove any instances of the current libs from the
        # run-time libraries.
        if "rt_lib_dirs" in cfg:
            new_lib_dirs = []
            for d in cfg["rt_lib_dirs"]:
                if d not in opts.get("lib_dir", []):
                    new_lib_dirs.append(d)
            cfg["rt_lib_dirs"] = new_lib_dirs
        if "rt_rpaths" in cfg:
            new_rpaths = []
            for d in cfg["rt_rpaths"]:
                if d not in opts.get("rpath", []):
                    new_rpaths.append(d)
            cfg["rt_rpaths"] = new_rpaths
        if "rt_libs" in cfg:
            new_libs = []
            for l in cfg["rt_libs"]:
                if l not in opts.get("lib", []):
                    new_libs.append(l)
            cfg["rt_libs"] = new_libs

        cfg.validate_modules(self)
        return True

    def test_options(self, cfg):
        opts = Option("mpi_show", "-show", type="bool")
        if not self.test_option(cfg, opts,
                                args=[("mpi_show", True)]):
            return False
        return True
