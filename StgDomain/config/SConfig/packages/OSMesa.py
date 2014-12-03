import os
import SConfig

class OSMesa(SConfig.Package):
    def __init__(self, scons_env, scons_opts, required=False):
        SConfig.Package.__init__(self, scons_env, scons_opts, required)
        self.libraries = [['GL', 'GLU']]
        self.have_define = 'HAVE_MESA'
