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
        self.help["location"] = "Specifies the command to use " + \
            "for compiling C code with the plan <name>."
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


    def parse_version_string(self, cfg):
        pass


#     def darwin_prog_args(self, cfg, env):
#         bak = env.get("rpath", None)
# 	env["rpath"] = []
#         cmd_args = self.conv_link_env(cfg, env)
# 	if bak is not None:
# 		env["rpath"] = bak
# 	return cmd_args


#     def darwin_shared_args(self, cfg, env):
#         cmd_args = []
#         if "install_name" in env["options"]:
#             cmd_args.append(("install_name", env["target"]))
#         cmd_args.extend([("libdir", d) for d in env.get("lib_dirs", []) \
#                              if d not in cfg["default_lib_dirs"]])
#         cmd_args.extend(env.get("lib_paths", []))
#         cmd_args.extend([("lib", d) for d in env.get("libs", [])])
#         return cmd_args


    def visit(self, cfg, visitor):
        visitor.ccompiler(self, cfg)
