import os
import SConfig

class libTIFF(SConfig.Package):
    def __init__(self, env, options):
        SConfig.Package.__init__(self, env, options)
        self.headers = [['tiff.h']]
        self.libraries = [['tiff']]
        self.have_define = 'HAVE_TIFF'
