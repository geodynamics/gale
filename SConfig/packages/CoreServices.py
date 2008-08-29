import os
import SConfig

class CoreServices(SConfig.Package):
    def __init__(self, scons_env, scons_opts, required=False, **kw):
        SConfig.Package.__init__(self, scons_env, scons_opts, required, **kw)
        self.frameworks = ['CoreServices']
        self.have_define = 'HAVE_CARBON'
