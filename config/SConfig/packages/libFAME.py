import os
import SConfig

class libFAME(SConfig.Package):
    def __init__(self, env, options):
        SConfig.Package.__init__(self, env, options)
        self.headers = [['fame.h']]
        self.libraries = [['fame']]
        self.have_define = 'HAVE_FAME'
