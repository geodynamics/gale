import os
import SConfig

class libJPEG(SConfig.Package):
    def __init__(self, env, options):
        SConfig.Package.__init__(self, env, options)
        self.headers = [['jpeglib.h']]
        self.libraries = [['jpeg']]
        self.have_define = 'HAVE_JPEG'
