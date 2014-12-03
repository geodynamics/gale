import os
import SConfig

class cmath(SConfig.Package):
    def __init__(self, scons_env, scons_opts, required=False):
        SConfig.Package.__init__(self, scons_env, scons_opts, required)
        self.dependency(SConfig.packages.CompilerFlags)
        self.libraries = [['m']]
