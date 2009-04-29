import os
import config
import config.utils as utils
from config.utils.options import Option
from CCompiler import CCompiler
from Linker import Linker

class gcc(CCompiler, Linker):

    def __init__(self, ctx):
        CCompiler.__init__(self, ctx)
        Linker.__init__(self, ctx)
        self.commands = ["gcc", "cc"]

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

    def apply_env(self, cfg, env):
        self.apply_compile_env(cfg, env)
        if cfg.is_valid(Linker):
            self.apply_link_env(cfg, env)

    def test_compile_options(self, cfg):
        lang = self.language
        src_empty = self.source_empty[lang]
        src_lib = self.source_lib[lang]

        # 32bit flag
        src_fn = self.ctx.make_temp_file(src_empty, ext=cfg["exts"][lang])
        out_fn = utils.path.replext(src_fn, ".o")
        opts = [Option("32bit", "-m32", type="bool"),
                Option("32bit", "-m32-bit", type="bool")]
        if not self.test_option(cfg, opts, products=out_fn,
                                args=[("compile", True), ("output", out_fn), ("32bit", True),
                                      src_fn]):
            return False

        # 64bit flag
        src_fn = self.ctx.make_temp_file(src_empty, ext=cfg["exts"][lang])
        out_fn = utils.path.replext(src_fn, ".o")
        opts = [Option("64bit", "-m64", type="bool"),
                Option("64bit", "-m64-bit", type="bool")]
        if not self.test_option(cfg, opts, products=out_fn,
                                args=[("compile", True), ("output", out_fn), ("64bit", True),
                                      src_fn]):
            return False

        # position independant code
        src_fn = self.ctx.make_temp_file(src_lib, ext=cfg["exts"][lang])
        out_fn = utils.path.replext(src_fn, ".o")
        opts = Option("pic", "-fPIC", type="bool")
        if not self.test_option(cfg, opts, products=out_fn,
                                args=[("compile", True), ("output", out_fn), ("pic", True),
                                      src_fn]):
            return False
        cfg["_pic_obj"] = out_fn

        # version
        opts = Option("version", "--version", type="bool")
        if not self.test_option(cfg, opts, products=out_fn,
                                args=[("version", True)]):
            return False

        # dump command-line info
        opts = Option("dump_cmds", "-###", type="bool")
        if not self.test_option(cfg, opts, err_okay=True,
                                args=[("dump_cmds", True), cfg["_simple_obj"]]):
            return False

        # flat namespace.
        out_fn = self.ctx.make_temp_file()
        opts = Option("flat_namespace", "-flat_namespace", type="bool")
        self.test_option(cfg, opts, "compile_cmd", products=out_fn,
                         args=[("output", out_fn),
                               ("flat_namespace", True),
                               src_fn])

        # single module
        out_fn = self.ctx.make_temp_file()
        opts = Option("single_module", "-single_module", type="bool")
        self.test_option(cfg, opts, "compile_cmd", products=out_fn,
                         args=[("output", out_fn),
                               ("single_module", True),
                               src_fn])

        # undefined symbols
        out_fn = self.ctx.make_temp_file()
        opts = Option("undefined", "-undefined", type="enum",
                      enum=["error", "warning", "suppress", "dynamic_lookup"])
        self.test_option(cfg, opts, "compile_cmd", products=out_fn,
                         args=[("output", out_fn),
                               ("undefined", "suppress"),
                               src_fn])

        return True

    def parse_version_string(self, cfg):
        res, stdout, stderr = self.execute(cfg, args=[("version", True)])
        if not res and not stderr:
            if stdout[:3] == "gcc":
                stdout = stdout.splitlines()[0]
                ver = stdout[stdout.rfind(' ') + 1:]
                cfg["version"] = ver.split('.')

    def get_link_rt_args(self, cfg, out, err):
        words = config.tools.Compiler.get_link_rt_args(self, cfg, out, err)

        # gcc puts quotes around everything, strip them.
        cmd_args = []
        for w in words:
            cmd_args.append(w[1:-1])
        return cmd_args
