import os
import SConfig

class dl(SConfig.Package):
    def __init__(self, env, options):
        SConfig.Package.__init__(self, env, options)
        self.headers = ['dlfcn.h']
        self.libraries = [['dl']]
