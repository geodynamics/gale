import os
import config
import config.utils as utils

class PETSc(config.Package):

    def __init__(self, ctx, **kw):
        config.Package.__init__(self, ctx, **kw)
        self.patterns = ["*petsc*", "*PETSc*", "PETSC"]
        self.inc_exts = ["petsc"]
        self.forced_arch = None

    def _setup_dependencies(self):
        config.Package._setup_dependencies(self)
        self.mpi = self.add_dependency(config.packages.MPI,
                                       required=True, combine=True)
        self.blas_lapack = self.add_dependency(config.packages.BlasLapack,
                                               required=True, combine=True)
        self.x11 = self.add_dependency(config.packages.X11,
                                       required=False, combine=True)


    def setup_libraries(self):
        incs = ["petsc.h", "petscvec.h", "petscmat.h",
                "petscksp.h", "petscsnes.h"]
        libs = ["petscsnes", "petscksp", "petscdm",
                "petscmat", "petscvec", "petsc"]
        self.add_library_set(incs, libs)
        self.add_auxilliary_libs("c", ["m"])


    def _gen_dir_sets(self, **kw):
        for ds in config.Package._gen_dir_sets(self, **kw):

            # If we have PETScExt as being configured, then we need to export the
            # base directory.
            if config.packages.PETScExt in self.ctx.modules:
                ds = (ds[0], ds[1] + [ds[0]], ds[2])

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


#     def gen_dirs(self):
#         utils.log("Searching for directories.")
#         utils.log.indent()
#         for base_dir in self.gen_base_dirs():
#             # Get the architecture.
#             archs = self.get_archs(base_dir)
#             # Insert an empty architecture at the beginning to
#             # find versions of petsc that don't need one.
#             archs.insert(0, None)
#             for a in archs:
#                 for set in self.gen_from_base_dir(base_dir, arch=a):
#                     yield set + [a]
#         utils.log.unindent()


#     def gen_lib_dir_sets(self, base_dir="", **kw):
#         arch = kw.get("arch", None)
#         for d in config.Package.gen_lib_dir_sets(self, base_dir):
#             # As far as I know PETSc should only need one of these.
#             if arch and len(d) == 1:
#                 # Uninstalled PETSc 3.
#                 yield [os.path.join(arch, d[0])]
#                 # PETSc 2.x
#                 yield [os.path.join(d[0], arch)]
#             else:
#                 yield d


#     def gen_inc_dir_sets(self, base_dir="", **kw):
#         arch = kw.get("arch", None)
#         for d in config.Package.gen_inc_dir_sets(self, base_dir):
#             yield d
#             if arch and len(d) == 1:
#                 # Uninstalled PETSc 3.
#                 yield [d[0], os.path.join(arch, d[0])]
#                 # PETSc 2.x
#                 yield [d[0], os.path.join("bmake", arch)]

#     def gen_dirs(self):
#         for dirs in config.Package.gen_dirs(self):
#             if not dirs[0] or len(dirs[1]) > 1 or len(dirs[2]) > 1:
#                 continue
#             archs = self.get_archs(dirs[0])
#             if archs:
#                 for a in archs:
#                     # PETSc 3 uninstalled
#                     yield [dirs[0], dirs[1],
#                            [os.path.join(dirs[0], a, "lib")],
#                            a]

#                     # PETSc 3 installed
#                     yield [dirs[0], dirs[1], dirs[2], a]

#                     # Other PETScs.
#                     yield [dirs[0],
#                            dirs[1] + [os.path.join(dirs[0], "bmake", a)],
#                            dirs[2] + [os.path.join(dirs[2][0], a)],
#                            a]
#             else:
#                 yield dirs + [None]

#     def candidate_from_dirs(self, dirs):
#         cfg = self.ctx.new_candidate(self, dirs[0], dirs[1], dirs[2])
#         cfg["arch"] = dirs[3]
#         return cfg

    def _scan_candidate(self, cfg):
        # Check if we have an architecture or not.
        arch = cfg.get("arch", None)
        if not arch:
            return

        # Add the bmake/arch include directory.
#         inc_dir = os.path.join('bmake', arch)
#         if os.path.exists(os.path.join(cfg.base_dir, inc_dir)):
#             self.add_inc_dirs(cfg, inc_dir)

        # Add the lib/arch library directory.
#         lib_dir = os.path.join('lib', cfg.arch)
#         if os.path.exists(os.path.join(cfg.base_dir, lib_dir)):
#             self.add_lib_dirs(cfg, lib_dir)

        # Parse extra libraries.
        extra_lib_dirs, extra_libs = self.get_extra_libraries(cfg)
#        cfg.add_lib_dirs(extra_lib_dirs)
#        cfg.add_aux_libs(extra_libs)

        # There seems to be some extra library paths we'll need that get stored
        # in SL_LINKER_LIBS.
#        extra_lib_dirs, extra_libs = self.get_sl_linker_libs(cfg)
#        cfg.add_lib_dirs(extra_lib_dirs)
#        cfg.add_libs(extra_libs)


    def _test(self, cfg):
        # Check for source code.
        if cfg["base_dir"]:
            p = os.path.join(cfg["base_dir"], "src", "mat", "impls", "aij", "mpi", "mpiaij.h")
            cfg["has_source"] = os.path.exists(p)
        else:
            cfg["has_source"] = False

        # Now business as usual.
        return config.Package._test(self, cfg)


    def make_run_command_line(self, cfg, fn, deps):
        mpirun = deps[0].get("mpirun", None)
        if mpirun is not None:
            return mpirun["make_command_line"](mpirun, args=[("nprocs", "1"),
                                                             os.path.join(".", fn)])
        return config.Package.make_run_command_line(self, cfg, fn, deps)


    def get_archs(self, base_dir):
        fnd_archs = []
        if base_dir:
            fn = os.path.join(base_dir, 'bmake', 'petscconf')
            if os.path.exists(fn):
                f = file(fn, 'r')
                a = f.readline().split('=')[1][:-1]
                f.close()
                if a not in fnd_archs:
                    fnd_archs.append(a)

            archs = os.listdir(base_dir)
            for a in archs:
                if os.path.isdir(os.path.join(base_dir, a)) \
                        and (a.find('opt') != -1 or a.find('debug') != -1):
                    fnd_archs.append(a)

        return fnd_archs


    def get_extra_libraries(self, cfg):
        """Read 'petscconf' and extract any additional dependencies/extra
        libraries we may need."""

        # Make sure the file exists before trying anything.
        petscconf = os.path.join(cfg["base_dir"], 'bmake', cfg["arch"], 'petscconf')
        if not os.path.exists(petscconf):
            return None, None

        # Read all the lines, which are of the form 'something = something else'.
        f = file(petscconf, 'r')
        line_dict = {}
        for line in f.readlines():
            sides = line.split('=')
            line_dict[sides[0].strip()] = sides[1].strip()
        f.close()

        # Try and locate any possible dependent cfgallations
        # PETSc knows about.
        name_map = {'MPI': self.mpi,
                    'BLASLAPACK': self.blas_lapack}
        for name_base, pkg in name_map.iteritems():
            name = name_base + '_INCLUDE'
            if name not in line_dict: continue
            string = utils.macro.subst(line_dict[name], line_dict).strip()
            for sub in string.split(' '):
                if sub[:2] == "-I":
                    base_dir = os.path.normpath(sub[2:])

                    # Try the base directory on it's own; sometimes
                    # the libraries will be placed there.
                    for cand in pkg.gen_cands_from_dir_set((base_dir, [base_dir], [base_dir])):
                        pkg.add_candidate(cand)

                    # Try combining with sub-directories.
                    base_dir = os.path.dirname(base_dir)
                    for ds in pkg.gen_from_base_dir(base_dir):
                        for cand in pkg.gen_cands_from_dir_set(ds):
                            pkg.add_candidate(cand)

            name = name_base + '_LIB'
            if name not in line_dict: continue
            string = utils.macro.subst(line_dict[name], line_dict).strip()
            for sub in string.split(' '):
                if sub[:2] == "-L":
                    base_dir = os.path.normpath(sub[2:])

                    # Try the base directory on it's own; sometimes
                    # the libraries will be placed there.
                    for cand in pkg.gen_cands_from_dir_set((base_dir, [base_dir], [base_dir])):
                        pkg.add_candidate(cand)

                    # Try combining with sub-directories.
                    base_dir = os.path.dirname(base_dir)
                    for ds in pkg.gen_from_base_dir(base_dir):
                        for cand in pkg.gen_cands_from_dir_set(ds):
                            pkg.add_candidate(cand)

        # Hunt down all the libraries and library paths we may need.
        names = ['PACKAGES_LIBS', 'SL_LINKER_LIBS']
        for name in names:
            if name not in line_dict: continue
            lib_string = line_dict['PACKAGES_LIBS']
            lib_string = utils.macro.subst(lib_string, line_dict)

            extra_lib_dirs = []
            extra_libs = []
            for string in lib_string.split(' '):
                if string[:2] == "-l":
                    extra_libs += [string[2:]]
                elif string[:2] == "-L":
                    extra_lib_dirs += [string[2:]]
        return (extra_lib_dirs, extra_libs)


    def get_sl_linker_libs(self, cfg):
        petscconf = os.path.join(cfg["base_dir"], 'bmake', cfg["arch"], 'petscconf')
        if not os.path.exists(petscconf): return ([], [])
        f = file(petscconf, 'r')
        line_dict = {}
        for line in f.readlines():
            sides = line.split('=')
            line_dict[sides[0].strip()] = sides[1].strip()
        f.close()
        if 'SL_LINKER_LIBS' not in line_dict: return ([], [])
        lib_string = line_dict['SL_LINKER_LIBS']
        lib_string = utils.macro.subst(lib_string, line_dict)

        extra_lib_dirs = []
        extra_libs = []
        for string in lib_string.split(' '):
            if string[:2] == "-l":
                extra_libs += [string[2:]]
            elif string[:2] == "-L":
                extra_lib_dirs += [string[2:]]
        return (extra_lib_dirs, extra_libs)


    source_code = {"c": """#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>
#include <petsc.h>
#include <petscvec.h>
#include <petscmat.h>
#include <petscksp.h>
int main( int argc, char** argv ) {
  Vec vec;
  Mat mat;
  KSP ksp;
  PetscInitialize( &argc, &argv, PETSC_NULL, PETSC_NULL );
  VecCreate( PETSC_COMM_WORLD, &vec );
  MatCreate( PETSC_COMM_WORLD, &mat );
  KSPCreate( PETSC_COMM_WORLD, &ksp );
  KSPDestroy( ksp );
  MatDestroy( mat );
  VecDestroy( vec );
  printf( "%d.%d.%d\\n", PETSC_VERSION_MAJOR, PETSC_VERSION_MINOR,
          PETSC_VERSION_SUBMINOR );
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
