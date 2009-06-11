import os
from Package import Package
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
        yield ('/usr/local', ['/usr/local/include/petsc'], [])
        yield ('/usr/local', ['/usr/local/include'], ['/usr/local/lib'])
        yield ('/usr/local', ['/usr/local/include/petsc'], ['/usr/local/lib'])

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
                    env.AppendUnique(CPPPATH=loc[1] + [os.path.dirname(petscconf)])
                    # Add arch to the library directory.
                    env.AppendUnique(LIBPATH=[os.path.join(loc[2][0], self.arch)])
                    env.AppendUnique(RPATH=[os.path.join(loc[2][0], self.arch)])
                    # Add additional libraries.
                    from distutils import sysconfig
                    vars = {}
                    sysconfig.parse_makefile(petscconf, vars)
                    flags = sysconfig.expand_makefile_vars(vars['PACKAGES_LIBS'], vars)
                    flag_dict = env.ParseFlags(flags)
                    if 'LIBS' in flag_dict:
                        extra_libs = flag_dict['LIBS']
                        del flag_dict['LIBS']
                    env.MergeFlags(flag_dict)
                    okay = True

        if not okay:
            env.AppendUnique(CPPPATH=loc[1])
            env.AppendUnique(LIBPATH=loc[2])
            env.AppendUnique(RPATH=loc[2])

        # Add the libraries.
        env.PrependUnique(LIBS=['petscsnes', 'petscksp', 'petscdm',
                                'petscmat', 'petscvec', 'petsc'])
        env.AppendUnique(LIBS=extra_libs)

        yield env

    def check(self, conf, env):
        return conf.CheckLibWithHeader(None,
                                       ['mpi.h', 'petsc.h', 'petscvec.h', 'petscmat.h',
                                        'petscksp.h', 'petscsnes.h'], 'c',
                                       autoadd=0)
