import os
import SConfig

class MPI(SConfig.Package):
    def __init__(self, env, options):
        SConfig.Package.__init__(self, env, options)
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

    def validate_location(self, location):
        for lib_dir in location[2]:
            shared_dir = os.path.join(lib_dir, 'shared')
            path = os.path.join(location[0], shared_dir)
            if os.path.exists(path):
                location[2] += [shared_dir]
        return [1, '', '']
