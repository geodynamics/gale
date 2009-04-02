import os
import config
from config import utils
from config.utils.options import Option

class MPI(config.Package, config.Tool):

    def __init__(self, ctx, **kw):
        config.Package.__init__(self, ctx, **kw)
        self.forced_mpirun = None


    def _setup_dependencies(self):
        # We have to create an instance of 'Tool' because we need
        # to give any mpirun candidates access to the 'execute' and
        # 'make_command_line' methods. We don't want to search
        # for tools though, so set the mode to 1 (disabled).
        config.Package._setup_dependencies(self)
        self.mpic = self.add_dependency(config.tools.MPICompiler, mode=0)
        self.tool = self.ctx.get_module(config.Tool)


    def _setup_members(self):
        if self.__module__ == "config.packages.MPI":
            self.add_member(config.packages.MPICH)
            self.add_member(config.packages.OpenMPI)


    def setup_trial(self, trial):
        config.Package.setup_trial(self, trial)
        trial["mpi_compilers"] = []


    def _scan_candidate(self, cand):
        # Check for compilers and mpirun.
        found = False
        for d in self.ctx.pkg_bin_dirs:
            bin_dir = os.path.join(cand["base_dir"], d)
            if not os.path.exists(bin_dir):
                continue

            path = os.path.join(bin_dir, "mpicc")
            if "mpicc" not in cand and os.path.exists(path):
                utils.log("Found %s" % path)
                cand["mpicc"] = self.ctx.new_candidate(
                    [config.tools.CCompiler, config.tools.MPICompiler],
                    path)
                cand.append_unique("mpi_compilers", cand, make_list=True)
                self.com.cc.add_candidate(cand["mpicc"])
                self.mpic.add_candidate(cand["mpicc"])
                found = True

            path = os.path.join(bin_dir, "mpif77")
            if "mpif77" not in cand and os.path.exists(path):
                utils.log("Found %s" % path)
                cand["mpif77"] = self.ctx.new_candidate(
                    [config.tools.F77Compiler, config.tools.MPICompiler],
                    path)
                cand.append_unique("mpi_compilers", cand, make_list=True)
                self.com.f77.add_candidate(cand["mpif77"])
                self.mpic.add_candidate(cand["mpif77"])
                found = True

            path = os.path.join(bin_dir, "mpif90")
            if "mpif90" not in cand and os.path.exists(path):
                utils.log("Found %s" % path)
                cand["mpif90"] = self.ctx.new_candidate(
                    [config.tools.F77Compiler, config.tools.MPICompiler],
                    path)
                cand.append_unique("mpi_compilers", cand, make_list=True)
                self.com.f77.add_candidate(cand["mpif90"])
                self.mpic.add_candidate(cand["mpif90"])
                found = True

            if "mpirun" not in cand:
                found = False
                if self.forced_mpirun is not None:
                    path = self.forced_mpirun
                    found = True
                else:
                    for cmd in ["mpirun", "mpiexec"]:
                        path = os.path.join(bin_dir, cmd)
                        if os.path.exists(path):
                            found = True
                            break
                if found:
                    utils.log("Found %s" % path)
                    cand["mpirun"] = self.ctx.new_candidate(config.Tool, path)
                    found = True
                    break
        return found


    def _test(self, cfg):
        if "mpirun" in cfg:
            if not self.test_mpirun(cfg["mpirun"])[0]:
                del cfg["mpirun"]
        res = config.Package._test(self, cfg)
        if res:
            cfg.validate_modules(MPI)
        return res


    def gen_compilers(self):
        for com in self.com.configs:
            if com not in self.mpic.configs or \
                    com in self._cfg["mpi_compilers"]:
                yield com


    def gen_linkers(self):
        for lnk in self.lnk.configs:
            try:
                if lnk not in self.mpic.configs or \
                       lnk in self._cfg["mpi_compilers"]:
                    yield lnk
            except:
                import pdb
                pdb.set_trace()


    def make_run_command_line(self, cfg, fn, deps):
        mpirun = cfg.get("mpirun", None)
        if mpirun is not None:
            return mpirun["make_command_line"](mpirun, args=[("nprocs", "1"),
                                                             os.path.join(".", fn)])
        return config.Package.make_run_command_line(self, cfg, fn, deps)


    def test_mpirun(self, mpirun):
        """Ensure we have a valid mpirun binary."""
        # It looks like we can test this without actually needing an
        # MPICH binary to do it, so grab any old binary from the list
        # of compilers.
        prog = None
        for c in self.com.configs:
            if "_simple_prog" in c:
                prog = os.path.join(".", self.com.configs[0]["_simple_prog"])
                break
        if not prog:
            # Can't test mpirun without a working program.
            return False, {}
        mpirun["make_command_line"] = self.tool.make_command_line
        mpirun["execute"] = self.tool.execute
        mpirun["success"] = self.success
        opts = [Option("nprocs", "-n"),
                Option("nprocs", "-np")]
        if not self.test_option(mpirun, opts,
                                args=[("nprocs", "1"), prog]):
            return False, {"err_code": self.result[0], "stdout": self.result[1],
                           "stderr": self.result[2]}
        return True, {}


    def parse_mpicc_show(self, cfg, trial):
        res, out, err = cfg.mpicc.execute(args=[("mpi_show", True)])
        args = cfg.mpicc.options.parse_option_string(out)
        d = cfg.mpicc.options.gather(args)

    def setup_options(self):
        config.Package.setup_options(self)


    source_code = {"c": """#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>
int main( int argc, char** argv ) {
  MPI_Comm comm;
  MPI_Init( &argc, &argv );
  MPI_Comm_dup( MPI_COMM_WORLD, &comm );
  printf( "%d.%d\\n", MPI_VERSION, MPI_SUBVERSION );
  MPI_Finalize();
  return EXIT_SUCCESS;
}
"""}


    def setup_options(self):
        config.Package.setup_options(self)
        self.options.add_option(
            utils.options.Option("mpirun", "--" + self.opt_name + "-run-command",
                                 pref_sep="="))

    def process_options(self):
        config.Package.process_options(self)
        mpirun = self.get_option("mpirun", None)
        if mpirun is not None:
            self.forced_mpirun = mpirun
