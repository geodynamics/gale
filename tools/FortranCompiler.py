import os
import config
import config.utils as utils
from Compiler import Compiler
from config.utils.options import Option
from Linker import Linker

class FortranCompiler(Compiler):

    source_empty={"f77": """      program add
      implicit none
      real a
      a = 2
      stop
      end
""", "f95": """program add
implicit none
real a
a = 2
stop
end
"""}
    source_lib={"f77": """      subroutine add
      implicit none
      real a
      a = 2
      return
      end
""", "f95": """subroutine add
implicit none
real a
a = 2
return
end
"""}
    suffixes={"f77": ".f", "f95": ".f95"}

    def __init__(self, ctx):
        Compiler.__init__(self, ctx)
        self.language = "f77"
        self.help["location"] = "Specifies the command to use " + \
            "for compiling fortran code with the plan <name>."
        self.setup_defaults()

    def _setup_members(self):
        self.ctx.get_module(Compiler).add_member(self)
        if self.__module__ == "config.tools.FortranCompiler":
            self.add_member(config.tools.fort77)
            self.add_member(config.tools.pgf77)
            self.add_member(config.tools.pgf95)

    def _test(self, cfg):
        if not self.test_core_options(cfg):
            return False
        if not self.test_compile(cfg):
            return False
        if not self.test_options(cfg):
            return False

        cfg.validate_modules([self, FortranCompiler, config.tools.Compiler])
        cfg["compile"] = self.compile
        cfg["compile_args"] = self.conv_compile_env
        cfg.append_unique("languages", "f77", make_list=True)
        cfg["success"] = self.success
        if self.test_link(cfg):
            cfg.validate_modules(Linker)
            cfg["link"] = self.link
            cfg["prog_args"] = self.prog_args
            cfg["shared_args"] = self.shared_args
        self.parse_version_string(cfg)
        self.parse_rt_libs(cfg)
        return True

    def test_core_options(self, cfg):
        # output filename
        src_fn = self.ctx.make_temp_file(self.source_empty[self.type],
                                         suffix=self.suffixes[self.type])
        out_fn = utils.path.remext(src_fn)
        opts = Option("output", "-o")
        if not self.test_option(cfg, opts, products=out_fn, err_okay=True,
                                args=[("output", out_fn), src_fn]):
            return False

        # compile only
        src_fn = self.ctx.make_temp_file(self.source_empty[self.type],
                                         suffix=self.suffixes[self.type])
        out_fn = utils.path.replext(src_fn, ".o")
        opts = Option("compile", "-c", type="bool")
        if not self.test_option(cfg, opts, products=out_fn, err_okay=True,
                                args=[("compile", True),
                                      ("output", out_fn),
                                      src_fn]):
            return False
        cfg["_simple_obj"] = out_fn

        # header directory
        src_fn = self.ctx.make_temp_file(self.source_empty[self.type],
                                         suffix=self.suffixes[self.type])
        out_fn = utils.path.replext(src_fn, ".o")
        opts = Option("incdir", "-I", sep=["", " "])
        if not self.test_option(cfg, opts, products=out_fn, err_okay=True,
                                args=[("compile", True),
                                      ("incdir", "."),
                                      ("output", out_fn),
                                      src_fn]):
            return False

        # CPP definition
        src_fn = self.ctx.make_temp_file(self.source_empty[self.type],
                                         suffix=self.suffixes[self.type])
        out_fn = utils.path.replext(src_fn, ".o")
        opts = Option("cpp_def", "-D", sep="")
        if not self.test_option(cfg, opts, products=out_fn, err_okay=True,
                                args=[("compile", True),
                                      ("cpp_def", "DUMMY"),
                                      ("output", out_fn),
                                      src_fn]):
            return False

        # debug only
        src_fn = self.ctx.make_temp_file(self.source_empty[self.type],
                                         suffix=self.suffixes[self.type])
        out_fn = utils.path.replext(src_fn, ".o")
        opts = Option("debug", "-g", type="bool")
        if not self.test_option(cfg, opts, products=out_fn, err_okay=True,
                                args=[("compile", True),
                                      ("debug", True),
                                      ("output", out_fn),
                                      src_fn]):
            return False

        return True

    def test_compile(self, cfg):
        src_fn = self.ctx.make_temp_file(self.source_empty[self.type],
                                         suffix=self.suffixes[self.type])
        out_fn = utils.path.replext(src_fn, ".o")
        res, out, err = self.compile(cfg, [("compile", True),
                                           ("output", out_fn),
                                           src_fn])
        return not res and os.path.exists(out_fn)

    def test_link(self, cfg):
        src_fn = self.ctx.make_temp_file(self.source_empty[self.type],
                                         suffix=self.suffixes[self.type])
        out_fn = utils.path.remext(src_fn)
        res, out, err = self.link(cfg, [("output", out_fn), src_fn])
        if res or not os.path.exists(out_fn):
            return False
        cfg["_simple_prog"] = out_fn
        return True

    def parse_version_string(self, cfg):
        pass

    def visit(self, cfg, visitor):
        visitor.fortrancompiler(self, cfg)
