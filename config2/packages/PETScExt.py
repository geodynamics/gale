import os
from config2 import Package
from PETSc import PETSc

class PETScExt(Package):

    def setup_dependencies(self):
        self.petsc = self.env.ConfigurePackage(PETSc, 'PETSc')

    def gen_locations(self):
        yield ('/usr/local', ['/usr/local/include'], ['/usr/local/lib'])

    def gen_envs(self, loc):
        env = self.env.Clone()

        # Must have the architecture as well.
        self.arch = self.petsc.arch

        if loc[1]:
            env.AppendUnique(CPPPATH=loc[1])
        if loc[2]:
            env.AppendUnique(LIBPATH=[os.path.join(loc[2][0], self.arch)])

        yield env

    def check(self, conf):
        return conf.CheckLibsWithHeader(['petscext_utils', 'petscext_snes', 'petscext_ksp',
                                         'petscext_pc', 'petscext_mat', 'petscext_vec',
                                         'petscext_helpers'],
                                        ['mpi.h', 'petsc.h', 'petscvec.h', 'petscmat.h',
                                         'petscksp.h', 'petscsnes.h',
                                         'petscext.h', 'petscext_vec.h', 'petscext_mat.h',
                                         'petscext_ksp.h', 'petscext_snes.h'], 'c')
