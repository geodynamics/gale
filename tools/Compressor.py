import os
import config
import config.utils as utils
from config.utils.options import Option

class Compressor(config.Tool):

    def __init__(self, ctx):
        config.Tool.__init__(self, ctx)


    def _setup_members(self):
        self.add_member(config.tools.tar)


    def setup_trial(self, cfg):
        config.Tool.setup_trial(self, cfg)
        cfg["compress"] = self.compress
        cfg["compress_cmd"] = "$command $output $source"
        cfg["extension"] = ""


    def _test(self, cfg):
        if not self.test_options(cfg):
            return False
        cfg.validate_modules([self, Compressor])
        return True


    def test_options(self, cfg):
        # 
        src_fn = self.ctx.make_temp_file(src_empty, ext=cfg["exts"][lang])
        opts = Option("compile", "-c", type="bool")
        if not self.test_option(cfg, opts,
                                args=[("compile", True), src_fn]):
            return False

        # output filename
        src_fn = self.ctx.make_temp_file(src_empty, ext=cfg["exts"][lang])
        out_fn = utils.path.replext(src_fn, ".o")
        opts = Option("output", "-o")
        if not self.test_option(cfg, opts, products=out_fn,
                                args=[("compile", True), ("output", out_fn),
                                      src_fn]):
            return False
        cfg["_simple_obj"] = out_fn

        # header directory
        src_fn = self.ctx.make_temp_file(src_empty, ext=cfg["exts"][lang])
        out_fn = utils.path.replext(src_fn, ".o")
        opts = Option("inc_dir", "-I", sep=["", " "])
        if not self.test_option(cfg, opts, products=out_fn,
                                args=[("compile", True), ("output", out_fn), ("inc_dir", "."),
                                      src_fn]):
            return False

        # CPP definition
        src_fn = self.ctx.make_temp_file(src_empty, ext=cfg["exts"][lang])
        out_fn = utils.path.replext(src_fn, ".o")
        opts = Option("pp_def", "-D", sep="")
        if not self.test_option(cfg, opts, products=out_fn,
                                args=[("compile", True), ("output", out_fn), ("pp_def", "DUMMY"),
                                      src_fn]):
            return False

        # debug only
        src_fn = self.ctx.make_temp_file(src_empty, ext=cfg["exts"][lang])
        out_fn = utils.path.replext(src_fn, ".o")
        opts = Option("debug", "-g", type="bool")
        if not self.test_option(cfg, opts, products=out_fn,
                                args=[("compile", True), ("output", out_fn), ("debug", True),
                                      src_fn]):
            return False

        return True


    def compress(self, cfg, lang="", env={}, **kw):
        cmd = "compress_cmd"
        return self.execute(cfg, cmd, env, **kw)
