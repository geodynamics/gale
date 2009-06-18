from config import Package

class MPI(Package):

    def gen_locations(self):
        yield ('', ['/usr/include/mpi'], [])
        yield ('', ['/usr/local/include'], ['/usr/local/lib'])
        yield ('', ['/usr/local/include/mpi'], ['/usr/local/lib'])

    def gen_envs(self, loc):
	# If we've been given an MPI compiler just try that.
	if self.env['CC'] in ['mpicc', 'mpicxx']:
            yield self.env.Clone()
            return

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
