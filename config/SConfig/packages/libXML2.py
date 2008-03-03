import os
import SConfig

class libXML2(SConfig.Package):
    def __init__(self, env, options):
        SConfig.Package.__init__(self, env, options)
        self.setup_search_defaults()
        self.setup_options()
        self.sub_dirs += [[[os.path.join('include', 'libxml2')], ['lib']]]
        self.headers = [os.path.join('libxml', 'parser.h')]
        self.libraries = [['xml2']]
