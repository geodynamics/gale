import os
import config
import config.utils as utils
from MPI import MPI

class OpenMPI(MPI):

    def __init__(self, ctx):
        MPI.__init__(self, ctx)
        self.patterns = ["*openmpi*", "*OpenMPI*"]
        self.inc_exts = ["mpi"]


    def _setup_dependencies(self):
        MPI._setup_dependencies(self)
        self.syml = self.add_dependency(config.tools.SymLister)


    def _setup_members(self):
        MPI._setup_members(self)
        self.ctx.get_module(MPI).add_member(self)


    def setup_libraries(self):
        self.add_library_set("mpi.h", ["mpi", "open-rte", "open-pal"])
        self.add_auxilliary_libs("c", ["nsl", "util"])


    def check_symbols(self, cfg, libs, lib_paths, shared):
        if not len(self.syml.configs):
            return False
        syml = self.syml.configs[0]
        args = [lib_paths["mpi"]]
        if shared and "dynamic" in syml["options"]:
            args.insert(0, ("dynamic", True))
        res, out, err = syml["list_symbols"](syml, args)
        if out.find("ompi") == -1:
            return False
        return True
