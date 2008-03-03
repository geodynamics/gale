import os
import SConfig

class SVNRevision(SConfig.Package):
    def __init__(self, env, options):
        SConfig.Package.__init__(self, env, options)
        self.checks = [self.extract_revision]
        self.define_name = 'VERSION'
        self.checkout_path = os.getcwd()

    def extract_revision(self):
        print 'Extracting subversion revision number ... ',
        svn_path = os.path.join(self.checkout_path, '.svn', 'entries')
        if not os.path.exists(svn_path):
            print '\n   Could not find .svn directory.'
            return [0, '', '']
        f = file(svn_path, 'r')
        f.readline()
        f.readline()
        f.readline()
        ver = self.env['ESCAPE']('"' + str(int(f.readline())) + '"')
        f.close()
        self.cpp_defines += [(self.define_name, ver)]
        print 'okay'
        return [1, '', '']
