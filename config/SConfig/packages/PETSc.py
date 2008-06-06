import os, re
import SConfig

class PETSc(SConfig.Package):
    def __init__(self, scons_env, scons_opts, required=False, **kw):
        SConfig.Package.__init__(self, scons_env, scons_opts, required, **kw)
        self.mpi = self.dependency(SConfig.packages.MPI)
        self.blas_lapack = self.dependency(SConfig.packages.BlasLapack)
        self.base_patterns = ['petsc*', 'PETSC*', 'PETSc*']
        self.header_sub_dir = ['petsc']
        self.headers = [['petsc.h',
                         'petscvec.h', 'petscmat.h',
                         'petscksp.h', 'petscsnes.h']]
        self.libraries = [['petscsnes', 'petscksp',
                           'petscmat', 'petscvec',
                           'petscdm', 'petsc',]]
        self.symbols = [(['PetscInitialize', 'PetscFinalize'], '')]
        self.symbol_calls = ['%s(&argc, &argv, NULL, NULL);', '%s();']

    def process_installation(self, inst):
        # Read the PETSc architecture.
        inst.arch = self.get_arch(inst.base_dir)
        if not inst.arch:
            return False

        # Add the bmake/arch include directory.
        hdr_dir = os.path.join('bmake', inst.arch)
        if not os.path.exists(os.path.join(inst.base_dir, hdr_dir)):
            return False # Can't continue without bmake include directory.
        inst.add_hdr_dirs(hdr_dir)

        # Add the lib/arch library directory.
        lib_dir = os.path.join('lib', inst.arch)
        if not os.path.exists(os.path.join(inst.base_dir, lib_dir)):
            return False # Must have correct library path.
        inst.add_lib_dirs(lib_dir)

        # Parse extra libraries.
        extra_lib_dirs, extra_libs = self.get_extra_libraries(inst)
        inst.add_lib_dirs(extra_lib_dirs)
        inst.add_extra_libs(extra_libs)

        # There seems to be some extra library paths we'll need that get stored
        # in SL_LINKER_LIBS.
#        extra_lib_dirs, extra_libs = self.get_sl_linker_libs(inst)
#        inst.add_lib_dirs(extra_lib_dirs)
#        inst.add_libs(extra_libs)

        # Everything's okay.
        return True

    def get_arch(self, base_dir):
        petscconf = os.path.join(base_dir, 'bmake',
                                 'petscconf')
        if not os.path.exists(petscconf):
            return None
        f = file(petscconf, 'r')
        arch = f.readline().split('=')[1][:-1]
        f.close()
        return arch

    def get_extra_libraries(self, inst):
        """Read 'petscconf' and extract any additional dependencies/extra
        libraries we may need."""

        # Make sure the file exists before trying anything.
        petscconf = os.path.join(inst.base_dir, 'bmake', inst.arch, 'petscconf')
        if not os.path.exists(petscconf):
            return

        # Read all the lines, which are of the form 'something = something else'.
        f = file(petscconf, 'r')
        line_dict = {}
        for line in f.readlines():
            sides = line.split('=')
            line_dict[sides[0].strip()] = sides[1].strip()
        f.close()

        # Try and locate any possible dependent installations
        # PETSc knows about.
        name_map = {'MPI': self.mpi,
                    'BLASLAPACK': self.blas_lapack}
        for name_base, pkg in name_map.iteritems():
            name = name_base + '_INCLUDE'
            if name not in line_dict: continue
            string = self.subst(line_dict[name], line_dict).strip()
            for sub in string.split(' '):
                if sub[:len(self.env['INCPREFIX'])] == self.env['INCPREFIX']:
                    base_dir = os.path.normpath(sub[len(self.env['INCPREFIX']):])

                    # Try the base directory on it's own; sometimes
                    # the libraries will be placed there.
                    pkg.add_candidate(SConfig.Installation(pkg, base_dir))

                    # Try combining with sub-directories.
                    base_dir = os.path.dirname(base_dir)
                    for inst in pkg.combine_base_dir(base_dir):
                        if pkg.process_installation(inst):
                            pkg.add_candidate(inst)
            name = name_base + '_LIB'
            if name not in line_dict: continue
            string = self.subst(line_dict[name], line_dict).strip()
            for sub in string.split(' '):
                if sub[:len(self.env['LIBDIRPREFIX'])] == self.env['LIBDIRPREFIX']:
                    base_dir = os.path.normpath(sub[len(self.env['LIBDIRPREFIX']):])

                    # Try the base directory on it's own; sometimes
                    # the libraries will be placed there.
                    pkg.add_candidate(SConfig.Installation(pkg, base_dir, [], ['']))

                    # Try combining with sub-directories.
                    base_dir = os.path.dirname(base_dir)
                    for inst in pkg.combine_base_dir(base_dir):
                        if pkg.process_installation(inst):
                            pkg.add_candidate(inst)

        # Hunt down all the libraries and library paths we may need.
        names = ['PACKAGES_LIBS', 'SL_LINKER_LIBS']
        for name in names:
            if name not in line_dict: continue
            lib_string = line_dict['PACKAGES_LIBS']
            lib_string = self.subst(lib_string, line_dict)

            extra_lib_dirs = []
            extra_libs = []
            for string in lib_string.split(' '):
                if string[:len(self.env['LIBLINKPREFIX'])] == self.env['LIBLINKPREFIX']:
                    extra_libs += [string[len(self.env['LIBLINKPREFIX']):]]
                elif string[:len(self.env['LIBDIRPREFIX'])] == self.env['LIBDIRPREFIX']:
                    extra_lib_dirs += [string[len(self.env['LIBDIRPREFIX']):]]
        return (extra_lib_dirs, extra_libs)

    def get_sl_linker_libs(self, inst):
        petscconf = os.path.join(inst.base_dir, 'bmake', inst.arch, 'petscconf')
        if not os.path.exists(petscconf): return ([], [])
        f = file(petscconf, 'r')
        line_dict = {}
        for line in f.readlines():
            sides = line.split('=')
            line_dict[sides[0].strip()] = sides[1].strip()
        f.close()
        if 'SL_LINKER_LIBS' not in line_dict: return ([], [])
        lib_string = line_dict['SL_LINKER_LIBS']
        lib_string = self.subst(lib_string, line_dict)

        extra_lib_dirs = []
        extra_libs = []
        for string in lib_string.split(' '):
            if string[:len(self.env['LIBLINKPREFIX'])] == self.env['LIBLINKPREFIX']:
                extra_libs += [string[len(self.env['LIBLINKPREFIX']):]]
            elif string[:len(self.env['LIBDIRPREFIX'])] == self.env['LIBDIRPREFIX']:
                extra_lib_dirs += [string[len(self.env['LIBDIRPREFIX']):]]
        return (extra_lib_dirs, extra_libs)

    def subst(self, line, line_dict):
        inp = [w.strip() for w in line.split()]
        out = []
        while len(inp):
            w = inp[0]
            inp = inp[1:]
            if self.is_macro(w):
                new_line = self.expand_macro(w, line_dict)
                new_words = [nw.strip() for nw in new_line.split()]
                inp = new_words + inp
            else:
                out += [w]
        return ' '.join(out)

    def expand_macro(self, macro, line_dict):
        if macro[:2] == '${' and macro[-1:] == '}':
            macro = macro[2:-1]
        elif macro[0] == '$':
            macro = macro[1:]
        if macro not in line_dict: return ''
        return line_dict[macro]

    def is_macro(self, word):
        if (word[:2] == '${' and word[-1:] == '}') or word[0] == '$':
            return True
        return False

    def get_check_headers_fail_reason(self, fail_logs):
        for log in fail_logs:
            if log.find('MPI_') != -1:
                return 'Selected MPI implementation incompatible.'
        return ''
