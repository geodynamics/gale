from config2 import Package

class MPI(Package):

    def gen_locations(self):
        yield ('', ['/usr/include/mpi'], [])
        yield ('', ['/usr/local/include'], ['/usr/local/lib'])
        yield ('', ['/usr/local/include/mpi'], ['/usr/local/lib'])

    def check(self, conf):
        # Not sure which extra libraries to check for, so try them all.
        extra_libs = [[], ['rt'], ['pthread', 'rt'],
                      ['dl'], ['dl', 'rt'], ['dl', 'pthread'],
                      ['dl', 'pthread', 'rt']]
        for libs in extra_libs:

            # Check for general MPI.
            if conf.CheckLibsWithHeader(['mpi'] + libs, 'mpi.h', 'c'):
                return True

            # Check for MPICH.
            if conf.CheckLibsWithHeader(['mpich'] + libs, 'mpi.h', 'c') or \
                    conf.CheckLibsWithHeader(['mpich', 'pmpich'] + libs, 'mpi.h', 'c'):
                return True

            # Check for OpenMPI.
            if conf.CheckLibsWithHeader(['mpi', 'open-rte', 'open-pal'], 'mpi.h', 'c') or \
                    conf.CheckLibsWithHeader(['mpi', 'open-rte', 'open-pal', 'nsl', 'util'],
                                             'mpi.h', 'c'):
                return True
