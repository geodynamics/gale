import os
import SConfig

class MPI(SConfig.Package):
    def __init__(self, scons_env, scons_opts, required=False, **kw):
        SConfig.Package.__init__(self, scons_env, scons_opts, required, **kw)
        self.dependency(SConfig.packages.CompilerFlags)
        self.base_patterns = ['mpich*', 'MPICH*']
        self.header_sub_dir = ['mpi']
        self.headers = [['mpi.h']]
        self.libraries = [['mpich'],
                          ['pmpich', 'mpich'],
                          ['mpich', 'rt'],
                          ['pmpich', 'mpich', 'rt'],
                          ['mpich', 'pvfs2'],
                          ['pmpich', 'mpich', 'pvfs2'],
                          ['mpich', 'rt', 'pvfs2'],
                          ['pmpich', 'mpich', 'rt', 'pvfs2'],
                          ['mpi'],
                          ['lam', 'mpi']]
        self.shared_libraries = ['mpich', 'pmpich', 'mpi', 'lam']
        self.extra_libraries = ['rt', 'pvfs2']
        self.symbols = [(['MPI_Init', 'MPI_Comm_dup', 'MPI_Finalize'], '')]
        self.symbol_setup = 'MPI_Comm comm_world;\n'
        self.symbol_calls = ['%s(&argc, &argv);',
                             '%s(MPI_COMM_WORLD, &comm_world);',
                             '%s();']
        self.init_code = 'MPI_Init(&argc, &argv);'
        self.fina_code = 'MPI_Finalize();'
        self.check_code = """int is_ready;
MPI_Initialized(&is_ready);
return is_ready;"""

    def process_installation(self, inst):
        SConfig.Package.process_installation(self, inst)

        # MPICH sometimes stores it's shared libraries in prefix/lib/shared.
        for lib_dir in inst.lib_dirs:
            shared_dir = os.path.join(lib_dir, 'shared')
            path = os.path.join(inst.base_dir, shared_dir)
            if os.path.exists(path):
                inst.add_lib_dirs(shared_dir, prepend=True)
