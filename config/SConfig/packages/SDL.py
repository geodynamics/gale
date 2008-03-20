import os
import SConfig

class SDL(SConfig.Package):
    def __init__(self, scons_env, scons_opts, required=False):
        SConfig.Package.__init__(self, scons_env, scons_opts, required)
        self.header_sub_dir = 'SDL'
        self.headers = [['SDL.h'],
                        ['SDL/SDL.h']] # For framework.
        self.libraries = [['SDL']]
        self.frameworks = [['SDL', 'Cocoa']]
