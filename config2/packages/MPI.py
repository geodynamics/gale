from Package import Package

class MPI(Package):

    def gen_locations(self):
        yield ('', ['/usr/include/mpi'], [])
        yield ('', ['/usr/local/include'], ['/usr/local/lib'])
        yield ('', ['/usr/local/include/mpi'], ['/usr/local/lib'])

    def gen_envs(self, loc):
        # First try an empty environment just in case we have
        # been given an MPI compiler.
        yield self.env.Clone()

        for env in Package.gen_envs(self, loc):

            # Not sure which extra libraries to check for, so try them all.
            extra_libs = [[], ['rt'], ['pthread', 'rt'],
                          ['dl'], ['dl', 'rt'], ['dl', 'pthread'],
                          ['dl', 'pthread', 'rt']]
            for libs in extra_libs:

                # Check for general MPI.
                if self.find_libraries(loc[2], 'mpi'):
                    lib_env = env.Clone()
                    lib_env.PrependUnique(LIBS=['mpi'] + libs)
                    yield lib_env

                # Check for MPICH.
                if self.find_libraries(loc[2], 'mpich'):
                    lib_env = env.Clone()
                    lib_env.PrependUnique(LIBS=['mpich'] + libs)
                    yield lib_env
                if self.find_libraries(loc[2], ['mpich', 'pmpich']):
                    lib_env = env.Clone()
                    lib_env.PrependUnique(LIBS=['pmpich', 'mpich'] + libs)
                    yield lib_env

                # Check for OpenMPI.
                if self.find_libraries(loc[2], 'mpi'):
                    lib_env = env.Clone()
                    lib_env.PrependUnique(LIBS=['mpi', 'open-rte', 'open-pal'] + libs)
                    yield lib_env
                if self.find_libraries(loc[2], 'mpi'):
                    lib_env = env.Clone()
                    lib_env.PrependUnique(LIBS=['mpi', 'open-rte', 'open-pal', 'nsl', 'util'] + libs)
                    yield lib_env

    def check(self, conf, env):
        return conf.CheckLibWithHeader(None, 'mpi.h', 'c', autoadd=0)
