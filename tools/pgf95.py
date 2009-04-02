import os
import config
import config.utils as utils
from config.utils.options import Option
from F77Compiler import F77Compiler
from Linker import Linker


class pgf95(F77Compiler, Linker):

    def __init__(self, ctx):
        F77Compiler.__init__(self, ctx)
        self.language = "f90"
        self.commands = ["pgf95"]

        self.source_empty["f90"] = """program add
implicit none
real a
a = 2
stop
end
"""
        self.source_lib["f90"] = """subroutine add
implicit none
real a
a = 2
return
end
"""


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
        opts = Option("dump_cmds", "-###", type="bool")
        if not self.test_option(cfg, opts, err_okay=True,
                                args=[("dump_cmds", True), cfg["_simple_obj"]]):
            return False

        # Make sure this isn't a PG compiler.
        src_fn = self.ctx.make_temp_file(src_lib, ext=cfg["exts"][lang])
        out_fn = utils.path.replext(src_fn, ".o")
        opts = Option("fast", "-fast", type="bool")
        if not self.test_option(cfg, opts, err_okay=True,
                                args=[("compile", True), ("output", out_fn), ("pic", True),
                                      ("fast", True),
                                      src_fn]):
            return False

        return True


    def test_option(self, cfg, opts, err_okay=True,
                    products=[], args=[]):
        return config.Tool.test_option(self, cfg, opts, err_okay, products, args)


    def parse_version_string(self, cfg):
        pass
