import os
import SConfig

class PETScExt(SConfig.Package):
    def __init__(self, env, options):
        SConfig.Package.__init__(self, env, options)
        self.pkg_petsc = self.dependency(SConfig.packages.PETSc)
        self.base_patterns = ['petscext*', 'PETSCEXT*', 'PETScExt*']
        self.header_sub_dir = 'petsc'
        self.headers = [['petscext.h',
                         'petscext_vec.h', 'petscext_mat.h',
                         'petscext_ksp.h', 'petscext_snes.h']]
        self.libraries = [['petscext_snes', 'petscext_ksp', 'petscext_pc',
                           'petscext_mat', 'petscext_vec',
                           'petscext_utils']]
        self.require_shared = True
        self.use_rpath = True
        self.have_define = 'HAVE_PETSCEXT'

    def validate_location(self, location):
        # Just use whatever architecture PETSc uses.
        arch = self.pkg_petsc.arch
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
