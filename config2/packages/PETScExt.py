import os
from Package import Package
from PETSc import PETSc

class PETScExt(Package):

    def setup_dependencies(self):
        self.petsc = self.add_dependency(PETSc)

    def gen_locations(self):
        yield ('/usr/local', ['/usr/local/include'], ['/usr/local/lib'])

    def gen_envs(self, loc):
        env = self.env.Clone()

        # Must have the architecture as well.
        self.arch = self.petsc.arch
        if self.arch is not None:
            if loc[2]:
                loc[2][0] = os.path.join(loc[2][0], self.arch)

        env.AppendUnique(CPPPATH=loc[1])
        env.AppendUnique(LIBPATH=loc[2])
        env.AppendUnique(RPATH=loc[2])

        # In addition, export the PETSc base directory.
        if self.petsc.location[0]:
            env.AppendUnique(CPPPATH=[self.petsc.location[0]])

        # Try each library set.
        libs = ['petscext_utils', 'petscext_snes', 'petscext_ksp',
                'petscext_pc', 'petscext_mat', 'petscext_vec',
                'petscext_helpers']
        if self.find_libraries(loc[2], libs):
            lib_env = env.Clone()
            lib_env.PrependUnique(LIBS=libs)
            yield lib_env
        libs = ['petscext_utils', 'petscext_snes', 'petscext_ksp',
                'petscext_pc', 'petscext_mat', 'petscext_vec']
        if self.find_libraries(loc[2], libs):
            lib_env = env.Clone()
            lib_env.PrependUnique(LIBS=libs)
            yield lib_env

    def check(self, conf, env):
        return conf.CheckLibWithHeader(None,
                                       ['mpi.h', 'petsc.h', 'petscvec.h', 'petscmat.h',
                                        'petscksp.h', 'petscsnes.h',
                                        'petscext.h', 'petscext_vec.h', 'petscext_mat.h',
                                        'petscext_ksp.h', 'petscext_snes.h'], 'c',
                                       autoadd=0)
