import os
import SConfig

class cmath(SConfig.Package):
    def __init__(self, env, options):
        SConfig.Package.__init__(self, env, options)
        self.libraries = [['m']]
