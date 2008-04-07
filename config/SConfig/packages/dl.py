import os
import SConfig

class dl(SConfig.Package):
    def __init__(self, scons_env, scons_opts, required=False):
        SConfig.Package.__init__(self, scons_env, scons_opts, required)
        self.dependency(SConfig.packages.CompilerFlags)
        self.headers = [['dlfcn.h']]
        self.libraries = [['dl']]

    def setup(self):
        SConfig.Node.setup(self)
