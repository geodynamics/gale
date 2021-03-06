import os
import SConfig

class libPNG(SConfig.Package):
    def __init__(self, scons_env, scons_opts, required=False):
        SConfig.Package.__init__(self, scons_env, scons_opts, required)
        self.headers = [['png.h']]
        self.libraries = [['png']]
        self.have_define = 'HAVE_PNG'
