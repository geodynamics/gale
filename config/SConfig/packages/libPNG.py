import os
import SConfig

class libPNG(SConfig.Package):
    def __init__(self, env, options):
        SConfig.Package.__init__(self, env, options)
        self.headers = [['png.h']]
        self.libraries = [['png']]
        self.have_define = 'HAVE_PNG'
