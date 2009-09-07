import os
from config import Package
from MPI import MPI

class PETSc(Package):

    def setup_dependencies(self):
        self.mpi = self.add_dependency(MPI)

    def setup_options(self):
        from SCons.Script.Main import AddOption
        Package.setup_options(self)
        AddOption('--petsc-arch', dest='petsc_arch', nargs=1, type='string',
                  action='store', help='PETSc architecture.')

    def gen_locations(self):
        yield ('/usr/local', ['/usr/local/include'], ['/usr/local/lib'])
        yield ('/usr/local', ['/usr/local/include'], ['/usr/local/lib64'])
        yield ('/usr/local', ['/usr/local/include/petsc'], ['/usr/local/lib'])
        yield ('/usr/local', ['/usr/local/include/petsc'], ['/usr/local/lib64'])

    def gen_envs(self, loc):
        env = self.env.Clone()

        # Must have the architecture as well.
        self.arch = self.get_option('petsc_arch')

        # Try to find PETSc information.
        extra_libs = []
        okay = False
        if loc[0]:
            bmake_dir = os.path.join(loc[0], 'bmake')
            # If we don't alrady have an arch, try read it.
            if not self.arch:
                petscconf = os.path.join(bmake_dir, 'petscconf')
                try:
                    inf = open(petscconf)
                    self.arch = inf.readline().split('=')[1].strip()
                except:
                    pass
            if self.arch is not None:
                petscconf = os.path.join(bmake_dir, self.arch, 'petscconf')
                # Does it exist?
                if os.path.exists(petscconf):
                    # Add the include directories.
                    loc[1].append(os.path.dirname(petscconf))
                    # Add arch to the library directory.
                    loc[2][0] = os.path.join(loc[2][0], self.arch)
                    # Add additional libraries.
                    from distutils import sysconfig
                    vars = {}
                    sysconfig.parse_makefile(petscconf, vars)
                    flags = sysconfig.expand_makefile_vars(vars['PACKAGES_LIBS'], vars)
                    if 'X11_INCLUDE' in vars:
                        flags += ' ' + sysconfig.expand_makefile_vars(str(vars['X11_INCLUDE']), vars)
                    if 'MPI_INCLUDE' in vars:
                        flags += ' ' + sysconfig.expand_makefile_vars(str(vars['MPI_INCLUDE']), vars)
                    if 'BLASLAPACK_INCLUDE' in vars:
                        flags += sysconfig.expand_makefile_vars(str(vars['BLASLAPACK_INCLUDE']), vars)
                    if 'PCC_LINKER_FLAGS' in vars:
                        flags += ' ' + sysconfig.expand_makefile_vars(str(vars['PCC_LINKER_FLAGS']), vars)
                    if 'PCC_FLAGS' in vars:
                        flags += ' ' + sysconfig.expand_makefile_vars(str(vars['PCC_FLAGS']), vars)
                    if 'PCC_LINKER_LIBS' in vars:
                        flags += ' ' + sysconfig.expand_makefile_vars(str(vars['PCC_LINKER_LIBS']), vars)
                    flag_dict = env.ParseFlags(flags)
                    if 'LIBS' in flag_dict:
                        extra_libs = flag_dict['LIBS']
                        del flag_dict['LIBS']
                    env.MergeFlags(flag_dict)
                    okay = True

        env.AppendUnique(CPPPATH=loc[1])
        env.AppendUnique(LIBPATH=loc[2])
        env.AppendUnique(RPATH=loc[2])

        # Add the libraries.
        libs = ['petscsnes', 'petscksp', 'petscdm',
                'petscmat', 'petscvec', 'petsc']
        if self.find_libraries(loc[2], libs):
            env.PrependUnique(LIBS=libs)
            env.AppendUnique(LIBS=extra_libs)
            yield env

    def check(self, conf, env):
        return conf.CheckLibWithHeader(None,
                                       ['mpi.h', 'petsc.h', 'petscvec.h', 'petscmat.h',
                                        'petscksp.h', 'petscsnes.h'], 'c',
                                       autoadd=0)
