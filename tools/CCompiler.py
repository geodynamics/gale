import os, platform
import config
import config.utils as utils
from Compiler import Compiler
from config.utils.options import Option
from Linker import Linker

class CCompiler(Compiler):

    def __init__(self, ctx, **kw):
        Compiler.__init__(self, ctx, **kw)
        self.language = "c"
        self.setup_defaults()
        self.source_empty["c"] = "int main(int argc, char* argv[]){return 0;}\n"
        self.source_lib["c"] = "int func() {return 0;}\n"

    def _setup_members(self):
        self.ctx.get_module(Compiler).add_member(self)
        if self.__module__ == "config.tools.CCompiler":
            self.add_member(config.tools.gcc)

    def setup_trial(self, cfg):
        Compiler.setup_trial(self, cfg)
        cfg["compile_cmd_c"] = cfg["compile_cmd"]

    def _test(self, cfg):
        if not Compiler._test(self, cfg):
            return False
        cfg.validate_modules(CCompiler)
        self.parse_version_string(cfg)
        return True

    def apply_compile_env(self, cfg, env):
        Compiler.apply_compile_env(self, cfg, env)
        env["ccom"] = cfg["command"]

    def parse_version_string(self, cfg):
        pass

    def setup_options(self):
        Compiler.setup_options(self)
        self.options.add_option(
            utils.options.Option("flags", "--cflags", pref_sep="=",
                                 help=self.help["flags"]))

    def process_options(self):
        Compiler.process_options(self)
        flags = self.get_option("flags", [])
        if flags:
            flags = [flags]
        self.compile_env.append_unique("dcomflags", make_list=True, *flags)

    def make_help_dict(self):
        Compiler.make_help_dict(self)
        self.help["location"] = "Specifies the command to use " + \
            "for compiling C code."
        self.help["flags"] = "Specify any additional flags for the C compiler."

    def visit(self, cfg, visitor):
        visitor.ccompiler(self, cfg)
