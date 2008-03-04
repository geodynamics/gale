import os
import SConfig

class MPI(SConfig.Package):
    def __init__(self, env, options):
        SConfig.Package.__init__(self, env, options)
        self.setup_search_defaults()
        self.setup_options()
        self.base_patterns = ['mpich*', 'MPICH*']
        self.header_sub_dir = 'mpi'
        self.headers = ['mpi.h']
        self.libraries = [['mpich'],
                          ['mpich', 'pmpich'],
                          ['mpich', 'rt'],
                          ['mpich', 'pmpich', 'rt'],
                          ['mpi']]
        self.require_shared = True
        self.use_rpath = True
