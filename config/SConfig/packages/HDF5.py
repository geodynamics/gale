import os
import SConfig

class HDF5(SConfig.Package):
    def __init__(self, env, options):
        SConfig.Package.__init__(self, env, options)
        self.headers = [['hdf5.h']]
        self.libraries = [['hdf5']]
        self.have_define = 'HAVE_HDF5'
