import os
import SConfig

class OSMesa(SConfig.Package):
    def __init__(self, scons_env, scons_opts, required=False, **kw):
        SConfig.Package.__init__(self, scons_env, scons_opts, required, **kw)
        self.header_sub_dir = ['GL']
        self.headers = ['osmesa.h']
        self.libraries = [['GL', 'GLU']]
