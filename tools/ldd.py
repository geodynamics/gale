import os
import config
import config.utils as utils
from config.utils.options import Option
from DLTool import DLTool

class ldd(DLTool):

    def __init__(self, ctx, **kw):
        DLTool.__init__(self, ctx, **kw)
        self.commands = ["ldd"]

    def _setup_members(self):
        self.ctx.get_module(DLTool).add_member(self)

    def _test(self, cfg):
        return self.test_list(cfg)

    def test_list(self, cfg):
        src_fn = self.com.configs[0]["_shared_lib"]
        res, out, err = cfg["list_libraries"](cfg, [src_fn])
        if res or err:
            return False
        syms = self.get_symbols(out)
        if not syms:
            return False
        return True

    def get_libraries(self, out):
        libs = {}
        for l in out.splitlines():
            words = l.split()
            if len(words) != 4:
                continue
            libs[words[0]] = words[2]
        return libs
