import os
import config
import config.utils as utils

class PETScExt(config.Package):

    def __init__(self, ctx):
        config.Package.__init__(self, ctx)
        self.patterns = ["*petscext*", "*PETScExt*", "PETSCEXT"]
        self.forced_arch = None


    def _setup_dependencies(self):
        config.Package._setup_dependencies(self)
        self.petsc = self.add_dependency(config.packages.PETSc,
                                         required=True, combine=True)


    def setup_libraries(self):
        incs = ["petscext.h", "petscext_vec.h", "petscext_mat.h",
                "petscext_ksp.h", "petscext_snes.h"]
        libs = ["petscext_utils", "petscext_snes", "petscext_ksp",
                "petscext_pc", "petscext_mat", "petscext_vec", "petscext_helpers"]
        self.add_library_set(incs, libs)
        self.add_library_set(incs, libs[:-1])


    def _gen_dir_sets(self, **kw):
        for ds in config.Package._gen_dir_sets(self, **kw):

            # Get the architecture.
            if self.forced_arch is not None:
                archs = [self.forced_arch]
            else:
                archs = self.get_archs(ds[0])

            # Insert an empty architecture at the beginning to
            # find versions of petsc that don't need one.
            archs.insert(0, None)
            for a in archs:

                if a:
                    yield (ds[0], # PETScExt 2.*
                           ds[1] + [os.path.join(ds[0], "bmake", a)],
                           [os.path.join(ds[0], "lib", a)],
                           a)
                    yield (ds[0], # PETScExt 3.* uninstalled
                           ds[1] + [os.path.join(ds[0], a, "include")],
                           [os.path.join(ds[0], a, "lib")],
                           a)

                else:
                    yield ds


    def finalise_candidate(self, dirs, info, cand):
        if len(dirs) == 4:
            cand["arch"] = dirs[3]


    def get_archs(self, base_dir):
        fnd_archs = []
        if base_dir:
            archs = os.listdir(base_dir)
            for a in archs:
                if os.path.isdir(os.path.join(base_dir, a)) \
                        and (a.find('opt') != -1 or a.find('debug') != -1):
                    fnd_archs.append(a)

            dir = os.path.join(base_dir, 'lib')
            if not os.path.exists(dir):
                return fnd_archs
            archs = os.listdir(dir)
            for a in archs:
                if os.path.isdir(os.path.join(dir, a)) and \
                        (a.find('opt') != -1 or a.find('debug') != -1) and \
                        a not in fnd_archs:
                    fnd_archs.append(a)
        return fnd_archs


    def gen_dependencies(self):
        # In order to build PETScExt we need the PETSc we're using
        # to have the source code installed.
        for deps in config.Package.gen_dependencies(self):
            if deps[0]["has_source"]:
                yield deps
            else:
                utils.log("Skipping dependencies as PETSc does not have source installed.")

    def apply_compile_dependency(self, dep, com, lang, env):
        res = config.Package.apply_compile_dependency(self, dep, com, lang, env)
        if res:
            if dep.is_valid(config.packages.PETSc):
                env.append_unique("inc_dirs", dep["base_dir"], make_list=True)
        return res

    def make_run_command_line(self, cfg, fn, deps):
        mpirun = deps[1].get("mpirun", None)
        if mpirun is not None:
            return mpirun["make_command_line"](mpirun, args=[("nprocs", "1"),
                                                os.path.join(".", fn)])
        return config.Package.make_run_command_line(self, cfg, fn, deps)


    source_code = {"c": """#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>
#include <petsc.h>
#include <petscvec.h>
#include <petscmat.h>
#include <petscksp.h>
#include <petscext_vec.h>
#include <petscext_mat.h>
#include <src/mat/impls/aij/mpi/mpiaij.h>
int main( int argc, char** argv ) {
  Vec vec;
  Mat mat;
  KSP ksp;
  PetscInitialize( &argc, &argv, PETSC_NULL, PETSC_NULL );
  PetscFinalize();
  return EXIT_SUCCESS;
}
"""}

    def setup_options(self):
        config.Package.setup_options(self)
        self.options.add_option(
            utils.options.Option("arch", "--" + self.opt_name + "-arch", pref_sep="="))

    def process_options(self):
        config.Package.process_options(self)
        arch = self.get_option("arch", None)
        if arch is not None:
            self.forced_arch = arch
