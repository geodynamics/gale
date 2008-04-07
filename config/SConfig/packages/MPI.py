import os
import SConfig

class MPI(SConfig.Package):
    def __init__(self, scons_env, scons_opts, required=False):
        SConfig.Package.__init__(self, scons_env, scons_opts, required)
        self.dependency(SConfig.packages.CompilerFlags)
        self.base_patterns = ['mpich*', 'MPICH*']
        self.header_sub_dir = 'mpi'
        self.headers = [['mpi.h']]
        self.libraries = [['mpich'],
                          ['mpich', 'pmpich'],
                          ['mpich', 'rt'],
                          ['mpich', 'pmpich', 'rt'],
                          ['mpi'],
                          ['lam', 'mpi']]
        self.shared_libraries = ['mpich', 'pmpich', 'mpi', 'lam']
        self.require_shared = True
        self.symbols = [(['MPI_Init', 'MPI_Finalize'], '')]
        self.symbol_calls = ['%s(&argc, &argv);', '%s();']

    def generate_locations(self):
        for loc in SConfig.Package.generate_locations(self):
            for lib_dir in loc[2]:
                shared_dir = os.path.join(lib_dir, 'shared')
                path = os.path.join(loc[0], shared_dir)
                if os.path.exists(path):
                    loc[2] = [shared_dir] + loc[2]
            yield loc
