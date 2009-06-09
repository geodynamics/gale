from config2 import Package

class MPI(Package):

    def gen_locations(self):
        yield ('', ['/usr/include/mpi'], [])
        yield ('', ['/usr/local/include'], ['/usr/local/lib'])
        yield ('', ['/usr/local/include/mpi'], ['/usr/local/lib'])

    def gen_envs(self, loc):
        for env in Package.gen_envs(self, loc):

            # Not sure which extra libraries to check for, so try them all.
            extra_libs = [[], ['rt'], ['pthread', 'rt'],
                          ['dl'], ['dl', 'rt'], ['dl', 'pthread'],
                          ['dl', 'pthread', 'rt']]
            for libs in extra_libs:

                # Check for general MPI.
                lib_env = env.Clone()
                lib_env.PrependUnique(LIBS=['mpi'] + libs)
                yield lib_env

                # Check for MPICH.
                lib_env = env.Clone()
                lib_env.PrependUnique(LIBS=['mpich'] + libs)
                yield lib_env
                lib_env = env.Clone()
                lib_env.PrependUnique(LIBS=['pmpich', 'mpich'] + libs)
                yield lib_env

                # Check for OpenMPI.
                lib_env = env.Clone()
                lib_env.PrependUnique(LIBS=['mpi', 'open-rte', 'open-pal'] + libs)
                yield lib_env
                lib_env = env.Clone()
                lib_env.PrependUnique(LIBS=['mpi', 'open-rte', 'open-pal', 'nsl', 'util'] + libs)
                yield lib_env

    def check(self, conf):
        return conf.CheckLibWithHeader(None, 'mpi.h', 'c')
