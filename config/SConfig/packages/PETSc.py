import os
import SConfig

class PETSc(SConfig.Package):
    def __init__(self, env, options):
        SConfig.Package.__init__(self, env, options)
        self.setup_search_defaults()
        self.setup_options()
        self.dependencies = [SConfig.packages.MPI]
        self.base_patterns = ['petsc*', 'PETSC*', 'PETSc*']
        self.sub_dirs += [[[os.path.join('include', 'petsc')], ['lib']]]
        self.headers = ['petsc.h',
                        'petscvec.h', 'petscmat.h',
                        'petscksp.h', 'petscsnes.h']
        self.libraries = [['petsc', 'petscdm',
                           'petscvec', 'petscmat',
                           'petscksp', 'petscsnes']]
        self.require_shared = True
        self.use_rpath = True
        self.have_define = 'HAVE_PETSC'

        # Other stuff.
        self.arch = ''

    def validate_location(self, location):
        # Must have a base directory.
        if not location[0]:
            return (1, '', '')

        # Get the architecture.
        arch = self.get_petsc_arch(location[0])
        if not arch:
            return (0, '', 'Could not read architecture from petscconf.')
        self.arch = arch

        # Add the bmake/arch include directory.
        hdr_dir = os.path.join('bmake', arch)
        if not os.path.exists(os.path.join(location[0], hdr_dir)):
            return [0, '', 'No bmake/<arch> directory.']
        if hdr_dir not in location[1]:
            location[1] += [hdr_dir]

        # Add the lib/arch library directory.
        if 'lib' in location[2]:
            location[2].remove('lib')
        lib_dir = os.path.join('lib', arch)
        if not os.path.exists(os.path.join(location[0], lib_dir)):
            return [0, '', 'No lib/<arch> directory.']
        if lib_dir not in location[2]:
            location[2] += [lib_dir]

        return [1, '', '']

    def get_headers_error_message(self, console):
        if console.find('MPI_') != -1:
            return 'Incompatible implementation of MPI.'
        return ''

    def get_petsc_arch(self, base_dir):
        petscconf = os.path.join(base_dir, 'bmake',
                                 'petscconf')
        if not os.path.exists(petscconf):
            return None
        f = file(petscconf, 'r')
        arch = f.readline().split('=')[1][:-1]
        f.close()
        return arch
