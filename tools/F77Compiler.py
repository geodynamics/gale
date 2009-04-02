import os
import config
import config.utils as utils
from Compiler import Compiler
from config.utils.options import Option
from Linker import Linker

class F77Compiler(Compiler):


    def __init__(self, ctx):
        Compiler.__init__(self, ctx)
        self.language = "f77"
        self.help["location"] = "Specifies the command to use " + \
            "for compiling fortran code with the plan <name>."
        self.setup_defaults()

        self.source_empty["f77"] = """      program add
      implicit none
      real a
      a = 2
      stop
      end
"""
        self.source_lib["f77"] = """      subroutine add
      implicit none
      real a
      a = 2
      return
      end
"""


    def _setup_members(self):
        self.ctx.get_module(Compiler).add_member(self)
        if self.__module__ == "config.tools.F77Compiler":
            self.add_member(config.tools.fort77)
            self.add_member(config.tools.gfortran)
            self.add_member(config.tools.pgf77)
            self.add_member(config.tools.pgf95)


    def setup_trial(self, cfg):
        Compiler.setup_trial(self, cfg)
        cfg["compile_cmd_f77"] = cfg["compile_cmd"]


    def _test(self, cfg):
        if not Compiler._test(self, cfg):
            return False
        cfg.validate_modules(F77Compiler)
        return True


    def parse_version_string(self, cfg):
        pass


    def visit(self, cfg, visitor):
        visitor.fortrancompiler(self, cfg)
