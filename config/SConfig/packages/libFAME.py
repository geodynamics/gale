import os
import SConfig

class libFAME(SConfig.Package):
    def __init__(self, scons_env, scons_opts, required=False):
        SConfig.Package.__init__(self, scons_env, scons_opts, required)
        self.headers = [['fame.h']]
        self.libraries = [['fame']]
        self.have_define = 'HAVE_FAME'
