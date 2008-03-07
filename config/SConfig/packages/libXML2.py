import os
import SConfig

class libXML2(SConfig.Package):
    def __init__(self, env, options):
        SConfig.Package.__init__(self, env, options)
        self.dependency(SConfig.packages.CompilerFlags)
        self.header_sub_dir = 'libxml2'
        self.headers = [[os.path.join('libxml', 'parser.h')]]
        self.libraries = [['xml2']]
