import os
import SConfig

class HDF5(SConfig.Package):
    def __init__(self, scons_env, scons_opts, required=False, **kw):
        SConfig.Package.__init__(self, scons_env, scons_opts, required, **kw)
        self.base_patterns = ['*hdf5*', '*HDF5*']
        self.headers = [['hdf5.h']]
        self.libraries = [['hdf5'],
                          ['hdf5', 'pthread'],
                          ['hdf5', 'pthread', 'z'],
                          ['hdf5', 'pthread', 'z', 'sz']]
        self.require_parallel = False

    def process_installation(self, inst):
        SConfig.Package.process_installation(self, inst)

        inst.parallel_support = False
        extra_lib_dirs = []
        extra_libs = []
        for lib_dir in inst.lib_dirs:
            set_file = os.path.join(inst.base_dir, lib_dir, 'libhdf5.settings')
            if os.path.exists(set_file):
                f = open(set_file, 'r')
                for line in f.readlines():
                    if line.find('Extra libraries') != -1:
                        dict = self.env.ParseFlags(line.split(':')[1])
                        extra_lib_dirs = dict.get('LIBPATH', [])
                        extra_libs = dict.get('LIBS', [])
                    if line.find('Parallel support') != -1:
                        psup = line.split(':')[1].strip()
                        inst.parallel_support = (psup == 'yes')
        disabled = """
        inst.add_lib_dirs(extra_lib_dirs)
        inst.add_extra_libs(extra_libs)"""

        # Do we need the math library?
        if 'm' in extra_libs:
            self.cmath = self.dependency(SConfig.packages.cmath)

        # Do we need to include the library for szip or pthread?
        if 'pthread' in extra_libs: inst.add_extra_libs('pthread')
        if 'z' in extra_libs: inst.add_extra_libs('z')
        if 'sz' in extra_libs: inst.add_extra_libs('sz')

        # If we have parallel support or we require parallel support,
        # add in MPI requirements.
        if inst.parallel_support or self.require_parallel:
            self.mpi = self.dependency(SConfig.packages.MPI)
            self.symbols = [(['H5Pset_dxpl_mpio', 'H5Pset_fapl_mpio'], '')]
            self.symbol_calls = ['%s(dxpl_props, H5FD_MPIO_COLLECTIVE);',
                                 '%s(fapl_props, MPI_COMM_WORLD, MPI_INFO_NULL);']
            self.symbol_setup = """hid_t dxpl_props, fapl_props;
MPI_Comm comm_world;

MPI_Init(&argc, &argv);
MPI_Comm_dup(MPI_COMM_WORLD, &comm_world);
dxpl_props = H5Pcreate(H5P_DATASET_XFER);
fapl_props = H5Pcreate(H5P_FILE_ACCESS);
"""
            self.symbol_teardown = """H5Pclose(dxpl_props);
H5Pclose(fapl_props);
MPI_Finalize();
"""

        return

    def get_check_symbols_fail_reason(self, fail_logs):
        for log in fail_logs:
            if log.find('_mpio\''):
                return 'Not a parallel HDF5 implementation.'
        return ''
