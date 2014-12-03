import os
import SConfig

class PETScExt(SConfig.Package):
    def __init__(self, scons_env, scons_opts, required=False):
        SConfig.Package.__init__(self, scons_env, scons_opts, required)
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

    def generate_locations(self):
        for loc in SConfig.Package.generate_locations(self):
            if not loc[0]:
                yield loc
                continue

            # Just use whatever architecture PETSc uses.
            arch = self.pkg_petsc.arch
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
