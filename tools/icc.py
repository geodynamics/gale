import os
import config
import config.utils as utils
from config.utils.options import Option
from CCompiler import CCompiler
from Linker import Linker

class icc(CCompiler, Linker):


    def __init__(self, ctx):
        CCompiler.__init__(self, ctx)
        self.commands = ["icc"]


    def setup_defaults(self):
        CCompiler.setup_defaults(self)
        Linker.setup_defaults(self)


    def _setup_dependencies(self):
        CCompiler._setup_dependencies(self)


    def _setup_members(self):
        self.ctx.get_module(CCompiler).add_member(self)
        self.ctx.get_module(config.tools.Linker).add_member(self)


    def setup_trial(self, cfg):
        CCompiler.setup_trial(self, cfg)
        Linker.setup_trial(self, cfg)


    def _test(self, cfg):
        if not CCompiler._test(self, cfg):
            return False
        if self.test_core_link_options(cfg) and self.test_link_options(cfg):
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
        opts = Option("pic", "-fpic", type="bool")
        if not self.test_option(cfg, opts, products=out_fn,
                                args=[("compile", True), ("output", out_fn), ("pic", True),
                                      src_fn]):
            return False
        cfg["_pic_obj"] = out_fn

        # version
        opts = Option("version", "-dumpversion", type="bool")
        if not self.test_option(cfg, opts, products=out_fn,
                                args=[("version", True)]):
            return False

        # dump command-line info
        opts = Option("dump_cmds", "-dryrun", type="bool")
        if not self.test_option(cfg, opts, err_okay=True,
                                args=[("dump_cmds", True), cfg["_simple_obj"]]):
            return False

        return True

    def parse_version_string(self, cfg):
        res, stdout, stderr = self.execute(cfg, args=[("version", True)])
        if not res and not stderr:
            if stdout[:3] == "icc":
                stdout = stdout.splitlines()[0]
                ver = stdout[stdout.rfind(' ') + 1:]
                cfg["version"] = ver.split('.')


    def get_link_rt_args(self, cfg, out, err):
        # icc spits out everything and newline's a lot of it.
        s = err.find("ld  \\")
        if s == -1:
            return []

        lines = err[s:].splitlines()
        cmd_args = []
        for l in lines:
            if l[-2:] == " \\":
                cmd_args.extend(l[:-2].split())
            else:
                cmd_args.extend(l.split())
                break

        return cmd_args
