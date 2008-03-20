import os
import SConfig

class PETSc(SConfig.Package):
    def __init__(self, scons_env, scons_opts, required=False):
        SConfig.Package.__init__(self, scons_env, scons_opts, required)
        self.dependency(SConfig.packages.MPI)
        self.base_patterns = ['petsc*', 'PETSC*', 'PETSc*']
        self.header_sub_dir = 'petsc'
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

            arch = self.get_petsc_arch(loc[0])
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

            yield loc

    def get_petsc_arch(self, base_dir):
        petscconf = os.path.join(base_dir, 'bmake',
                                 'petscconf')
        if not os.path.exists(petscconf):
            return None
        f = file(petscconf, 'r')
        arch = f.readline().split('=')[1][:-1]
        f.close()
        return arch

    def get_check_headers_fail_reason(self, fail_logs):
        for log in fail_logs:
            if log.find('MPI') != -1:
                return 'Selected MPI implementation incompatible.'
        return ''
