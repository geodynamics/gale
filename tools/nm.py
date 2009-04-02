import os
import config
import config.utils as utils
from config.utils.options import Option
from SymLister import SymLister

class nm(SymLister):

    def __init__(self, ctx):
        SymLister.__init__(self, ctx)
        self.commands = ["nm"]

    def _setup_members(self):
        self.ctx.get_module(SymLister).add_member(self)

    def _test(self, cfg):
        if self.test_list(cfg):
            cfg.validate_modules([self, SymLister])
            cfg["list_symbols"] = self.list_symbols
            cfg["get_symbols"] = self.get_symbols
            self.test_dynamic(cfg)
            return True
        return False

    def test_list(self, cfg):
        src_fn = self.comp.configs[0]["_simple_obj"]
        args = [src_fn]
        res, out, err = self.execute(cfg, args)
        if res or err:
            return False
        syms = self.get_symbols(out)
        if not syms:
            return False
        return True

    def test_dynamic(self, cfg):
        src_fn = self.comp.configs[0]["_shared_lib"]
        args = [src_fn]
        opts = Option("dynamic", "-D", type="bool")
        if not self.test_option(cfg, opts,
                                args=[("dynamic", True), src_fn]):
            return False
        return True

    def get_symbols(self, out):
        syms = {}
        for l in out.splitlines():
            words = l.split()
            if len(words) < 2:
                continue
            if words[-2] in ["T", "s", "S"]:
                if "defined" not in syms:
                    syms["defined"] = []
                syms["defined"].append(words[-1])
            elif words[-2] in ["U", "u"]:
                if "undefined" not in syms:
                    syms["undefined"] = []
                syms["undefined"].append(words[-1])
        return syms
