import os
import SConfig

class dl(SConfig.Package):
    def __init__(self, env, options):
        SConfig.Package.__init__(self, env, options)
        self.dependency(SConfig.packages.CompilerFlags)
        self.headers = [['dlfcn.h']]
        self.libraries = [['dl']]
