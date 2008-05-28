import os
import SConfig

class HDF5(SConfig.Package):
    def __init__(self, scons_env, scons_opts, required=False, **kw):
        SConfig.Package.__init__(self, scons_env, scons_opts, required, **kw)
        self.require_parallel = False
        self.headers = [["hdf5.h"]]
        self.libraries = [["hdf5"]]

    def setup(self):
        SConfig.Package.setup(self)
        if self.require_parallel:
            self.dependency(SConfig.packages.MPI)
            self.symbols = [(["H5Pset_dxpl_mpio", "H5Pset_fapl_mpio"], "")]
            self.symbol_calls = ["%s( dxpl_props, H5FD_MPIO_COLLECTIVE );",
                                 "%s( fapl_props, MPI_COMM_WORLD, MPI_INFO_NULL );"]
            self.symbol_setup = """hid_t dxpl_props, fapl_props;
MPI_Init( &argc, &argv );
dxpl_props = H5Pcreate( H5P_DATASET_XFER );
fapl_props = H5Pcreate( H5P_FILE_ACCESS );
"""
            self.symbol_teardown = """H5Pclose( dxpl_props );
H5Pclose( fapl_props );
MPI_Finalize();
"""

    def generate_locations(self):
        for loc in SConfig.Package.generate_locations(self):
            extra_lib_dirs=[]
            extra_libs = []
            for lib_dir in loc[2]:
                set_file = os.path.join(loc[0], lib_dir, "libhdf5.settings")
                if os.path.exists(set_file):
                    f = open(set_file, "r")
                    for line in f.readlines():
                        if line.find("Extra libraries") != -1:
                            dict = self.env.ParseFlags(line.split(":")[1])
                            extra_libs = dict.get('LIBS', [])
                            extra_lib_dirs = dict.get('LIBPATH', [])
            old_libs = self.extra_libraries
            self.extra_libraries += extra_libs
            for lib_dir in extra_lib_dirs:
                if lib_dir not in loc[2]:
                    loc[2] += [lib_dir]
            yield loc
            self.extra_libraries = old_libs

    def get_check_symbols_fail_reason(self, fail_logs):
        for log in fail_logs:
            if log.find("_mpio'"):
                return "Not a parallel HDF5 implementation."
        return ''
