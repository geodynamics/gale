import os
import config
import config.utils as utils
from MPI import MPI

class MPICH(MPI):

    def __init__(self, ctx):
        MPI.__init__(self, ctx)
        self.patterns = ["*mpich*", "*MPICH*"]
        self.inc_exts = ["mpi"]
        self.no_mpd = False


    def _setup_dependencies(self):
        MPI._setup_dependencies(self)
        self.syml = self.add_dependency(config.tools.SymLister)


    def _setup_members(self):
        MPI._setup_members(self)
        self.ctx.get_module(MPI).add_member(self)


    def setup_libraries(self):
        self.add_library_set("mpi.h", "mpich")
        self.add_library_set("mpi.h", ["pmpich", "mpich"])
        self.add_library_set("mpi.h", "mpi")

        self.add_auxilliary_libs("c", ["rt"])
        self.add_auxilliary_libs("c", ["pthread"])
        self.add_auxilliary_libs("c", ["pthread", "rt"])
        self.add_auxilliary_libs("c", ["dl", "rt"])
        self.add_auxilliary_libs("c", ["dl", "pthread"])
        self.add_auxilliary_libs("c", ["dl", "pthread", "rt"])


    def check_symbols(self, cfg, libs, lib_paths, shared):
        if not len(self.syml.configs):
            raise "Need symbol lister for MPICH."
        syml = self.syml.configs[0]
        if "mpich" in lib_paths:
            args = [lib_paths["mpich"]]
        else:
            args = [lib_paths["mpi"]]
        if shared and "dynamic" in syml["options"]:
            args.insert(0, ("dynamic", True))
        res, out, err = syml["list_symbols"](syml, args)
        if out.find("ompi") != -1:
            return False
        return True


    def test_mpirun(self, mpirun):
        res = config.packages.MPI.test_mpirun(self, mpirun)
        if not self.no_mpd and "stderr" in res[1] and \
                res[1]["stderr"].find("no mpd is running") != -1:
            self.no_mpd = True
            self.messages.append("Warning: Couldn't find a valid 'mpirun' or 'mpiexec', " + \
                                     "seemingly because there is no 'mpd' running on this " + \
                                     "machine. If no MPI candidates were found please " + \
                                     "start 'mpd' and rerun the configuration.")
        return res
