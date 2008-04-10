import os, re
import SConfig

class PETSc(SConfig.Package):
    def __init__(self, scons_env, scons_opts, required=False, **kw):
        SConfig.Package.__init__(self, scons_env, scons_opts, required, **kw)
        self.dependency(SConfig.packages.MPI)
        self.base_patterns = ['petsc*', 'PETSC*', 'PETSc*']
        self.header_sub_dir = ['petsc']
        self.headers = [['petsc.h',
                         'petscvec.h', 'petscmat.h',
                         'petscksp.h', 'petscsnes.h']]
        self.libraries = [['petscsnes', 'petscksp',
                           'petscmat', 'petscvec',
                           'petscdm', 'petsc',]]
        self.require_shared = True
        self.symbols = [(['PetscInitialize', 'PetscFinalize'], '')]
        self.symbol_calls = ['%s(&argc, &argv, NULL, NULL);', '%s();']

        # Will be set after configuration.
        self.arch = ''

    def generate_locations(self):
        for loc in SConfig.Package.generate_locations(self):
            if not loc[0]:
                yield loc
                continue

            arch = self.get_arch(loc[0])
            if not arch:
                yield loc
                continue
            self.arch = arch

            # Add the bmake/arch include directory.
            hdr_dir = os.path.join('bmake', arch)
            if not os.path.exists(os.path.join(loc[0], hdr_dir)):
                continue
            if hdr_dir not in loc[1]:
                loc[1] += [hdr_dir]

            # Add the lib/arch library directory.
            if 'lib' in loc[2]:
                loc[2].remove('lib')
            lib_dir = os.path.join('lib', arch)
            if not os.path.exists(os.path.join(loc[0], lib_dir)):
                continue
            if lib_dir not in loc[2]:
                loc[2] += [lib_dir]

            # Parse extra libraries.
            extra_lib_dirs, extra_libs = self.get_extra_libraries(loc[0])
            if extra_lib_dirs: loc[2] += extra_lib_dirs
            if extra_libs: self.extra_libraries += extra_libs

            yield loc

    def get_arch(self, base_dir):
        petscconf = os.path.join(base_dir, 'bmake',
                                 'petscconf')
        if not os.path.exists(petscconf):
            return None
        f = file(petscconf, 'r')
        arch = f.readline().split('=')[1][:-1]
        f.close()
        return arch

    def get_extra_libraries(self, base_dir):
        petscconf = os.path.join(base_dir, 'bmake', self.arch, 'petscconf')
        if not os.path.exists(petscconf): return ([], [])
        f = file(petscconf, 'r')
        line_dict = {}
        for line in f.readlines():
            sides = line.split('=')
            line_dict[sides[0].strip()] = sides[1].strip()
        f.close()
        if 'PACKAGES_LIBS' not in line_dict: return ([], [])
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
