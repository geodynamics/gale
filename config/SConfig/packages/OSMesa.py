import os
import SConfig

class OSMesa(SConfig.Package):
    def __init__(self, env, options):
        SConfig.Package.__init__(self, env, options)
        self.libraries = [['GL', 'GLU']]
        self.have_define = 'HAVE_MESA'
