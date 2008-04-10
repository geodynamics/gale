import os
import SConfig

class libFAME(SConfig.Package):
    def __init__(self, scons_env, scons_opts, required=False, **kw):
        SConfig.Package.__init__(self, scons_env, scons_opts, required, **kw)
        self.headers = [['fame.h']]
        self.libraries = [['fame']]
        self.have_define = 'HAVE_FAME'
