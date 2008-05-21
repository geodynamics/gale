import os
import SConfig

class szip(SConfig.Package):
    def __init__(self, scons_env, scons_opts, required=False, **kw):
        SConfig.Package.__init__(self, scons_env, scons_opts, required, **kw)
        self.headers = [['szlib.h']]
        self.libraries = [['sz']]
