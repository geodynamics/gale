import os
import SConfig

class HDF5(SConfig.Package):
    def __init__(self, scons_env, scons_opts, required=False):
        SConfig.Package.__init__(self, scons_env, scons_opts, required)
        self.headers = [['hdf5.h']]
        self.libraries = [['hdf5']]
