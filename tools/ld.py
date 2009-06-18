import os
import config
import config.utils as utils
from config.utils.options import Option
from Linker import Linker

class ld(Linker):


    def __init__(self, ctx):
        Linker.__init__(self, ctx)
        self.commands = ["ld"]


    def _setup_members(self):
        self.ctx.get_module(config.tools.Linker).add_member(self)


    def setup_trial(self, cfg):
        Linker.setup_trial(self, cfg)

        # ld requires run-time symbols for 'main', or 'start' or whatever the
        # hell the compiler used to compile the objects being linked needs.
#         cfg["_rt_lib_dirs"] = "$string(lib_dir,$rt_lib_dirs)"
#         cfg["_rt_rpaths"] = "$string(rpath,$rt_rpaths)"
#         cfg["_rt_libs"] = "$string(lib,$rt_libs)"
#         cfg["_rt_begin"] = "$rt_begin"
#         cfg["_rt_end"] = "$rt_end"

        # Also needs this flag for creating the program header.
        cfg.append_unique("lnkprogflags", "-m", "elf_i386", "--eh-frame-hdr",
                          "-z", "relro", "--hash-style=both")

        # Looks like the dynamic linker will depend on the operating system.
        cfg.append_unique("lnkprogflags", "-dynamic-linker", "/lib/ld-linux.so.2")
