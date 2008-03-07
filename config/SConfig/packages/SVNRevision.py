import os
import SConfig

class SVNRevision(SConfig.Package):
    def __init__(self, env, options):
        SConfig.Package.__init__(self, env, options)
        self.checks = [self.extract_revision]
        self.define_name = 'VERSION'
        self.checkout_path = os.getcwd()

    def setup_options(self):
        pass

    def extract_revision(self):
        svn_path = os.path.join(self.checkout_path, '.svn', 'entries')
        if not os.path.exists(svn_path):
            return [0, '', 'Could not find .svn directory']
        f = file(svn_path, 'r')
        f.readline()
        f.readline()
        f.readline()
        ver = self.env['ESCAPE']('"' + str(int(f.readline())) + '"')
        f.close()
        self.cpp_defines += [(self.define_name, ver)]
        return [1, '', '']
