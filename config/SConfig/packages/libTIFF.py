import os
import SConfig

class libTIFF(SConfig.Package):
    def __init__(self, scons_env, scons_opts, required=False):
        SConfig.Package.__init__(self, scons_env, scons_opts, required)
        self.headers = [['tiff.h']]
        self.libraries = [['tiff']]
        self.have_define = 'HAVE_TIFF'
