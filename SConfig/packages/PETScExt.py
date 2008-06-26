import os
import SConfig

class PETScExt(SConfig.Package):
    def __init__(self, scons_env, scons_opts, required=False):
        SConfig.Package.__init__(self, scons_env, scons_opts, required)
        self.mpi = self.dependency(SConfig.packages.MPI)
        self.petsc = self.dependency(SConfig.packages.PETSc)
        self.base_patterns = ['petscext*', 'PETSCEXT*', 'PETScExt*']
        self.header_sub_dir = ['petsc']
        self.headers = [['petscext.h',
                         'petscext_vec.h', 'petscext_mat.h',
                         'petscext_ksp.h', 'petscext_snes.h']]
        self.libraries = [['petscext_utils',
                           'petscext_snes', 'petscext_ksp', 'petscext_pc',
                           'petscext_mat', 'petscext_vec']]

    def process_installation(self, inst):
        # Have to get the architecture using this garbage...
        archs = os.listdir(os.path.join(inst.base_dir, 'lib'))
        for arch in archs:
            if arch[0] != '.':
                inst.arch = arch

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

        return True
