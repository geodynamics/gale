import os
import config
import config.utils as utils
from config.utils.options import Option
from F77Compiler import F77Compiler
from Linker import Linker

class fort77(F77Compiler, Linker):


    def __init__(self, ctx):
        F77Compiler.__init__(self, ctx)
        Linker.__init__(self, ctx)
        self.commands = ["f77", "fort77"]


    def setup_defaults(self):
        F77Compiler.setup_defaults(self)
        Linker.setup_defaults(self)


    def _setup_dependencies(self):
        F77Compiler._setup_dependencies(self)


    def _setup_members(self):
        self.ctx.get_module(F77Compiler).add_member(self)
        self.ctx.get_module(config.tools.Linker).add_member(self)


    def setup_trial(self, cfg):
        F77Compiler.setup_trial(self, cfg)
        Linker.setup_trial(self, cfg)


    def setup_trial(self, cfg):
        F77Compiler.setup_trial(self, cfg)
        Linker.setup_trial(self, cfg)
        cfg["compile_success"] = self.compile_success


    def _test(self, cfg):
        if not F77Compiler._test(self, cfg):
            return False
        if self.test_link_options(cfg):
            cfg.validate_modules(Linker)
        self.parse_rt_libs(cfg)
        return True


    def test_compile_options(self, cfg):
        lang = self.language
        src_empty = self.source_empty[lang]
        src_lib = self.source_lib[lang]

        # position independant code
        src_fn = self.ctx.make_temp_file(src_lib, ext=cfg["exts"][lang])
        out_fn = utils.path.replext(src_fn, ".o")
        opts = Option("pic", "-fPIC", type="bool")
        if not self.test_option(cfg, opts, products=out_fn,
                                args=[("compile", True), ("output", out_fn), ("pic", True),
                                      src_fn]):
            return False
        cfg["_pic_obj"] = out_fn

        # dump command-line info
        opts = Option("dump_cmds", "-Wl,-###", type="bool")
        if not self.test_option(cfg, opts, err_okay=True,
                                args=[("dump_cmds", True), cfg["_simple_obj"]]):
            return False

        return True

    def test_link_options(self, cfg):
        if not Linker.test_link_options(self, cfg):
            return False

        # rpath
        out_fn = utils.path.remext(cfg["_simple_obj"])
        opts = Option("rpath", "-Wl,-Xlinker -Wl,-rpath -Wl,-Xlinker -Wl,", sep=[''])
        if not self.test_option(cfg, opts, products=out_fn,
                                args=[("output", out_fn), ("rpath", os.getcwd())] + \
                                    cfg.get("_static_start_objs", []) + \
                                    [cfg["_simple_obj"]] + cfg.get("_static_end_objs", []) + \
                                    cfg.get("_static_rt_args", [])):
            return False
        return True

    def test_option(self, cfg, opts, err_okay=True,
                    products=[], args=[]):
        return config.Tool.test_option(self, cfg, opts, err_okay, products, args)

    def get_link_rt_args(self, cfg, out, err):
        words = config.tools.Compiler.get_link_rt_args(self, cfg, out, err)

        # fort77 uses gcc, which puts quotes around everything, strip them.
        cmd_args = []
        for w in words:
            cmd_args.append(w[1:-1])
        return cmd_args

    def compile_success(self, res, out, err, err_okay=False, products=[]):
        return config.Tool.success(self, res, out, err, True, products)

    def parse_version_string(self, cfg):
        pass
