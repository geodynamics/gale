import os
import SConfig

class SDL(SConfig.Package):
    def __init__(self, env, options):
        SConfig.Package.__init__(self, env, options)
        self.header_sub_dir = 'SDL'
        self.headers = [['SDL.h'],
                        ['SDL/SDL.h']] # For framework.
        self.libraries = [['SDL']]
        self.frameworks = [['SDL', 'Cocoa']]
        self.have_define = 'HAVE_SDL'
