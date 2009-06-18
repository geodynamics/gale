import os
import config
import config.utils as utils
from config.utils.options import Option
from Compressor import Compressor


class tar(Compressor):

    def __init__(self, ctx):
        Compressor.__init__(self, ctx)
        self.commands = ["tar"]


    def _setup_members(self):
        self.ctx.get_module(Compressor).add_member(self)


    def test_options(self, cfg):
        garbage = "garbage" * 20

        # compress using gzip
        src_fn = self.ctx.make_temp_file(garbage)
        out_fn = src_fn + ".tgz"
        opts = Option("compress_output", "czf", type="str", sep=" ")
        if self.test_option(cfg, opts, products=out_fn,
                            args=[("compress_output", out_fn), src_fn]):
            cfg["extension"] = ".tgz"
            cfg["output"] = "$string(compress_output,$target)"
            return True

        return False
