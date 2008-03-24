import os
import SConfig

class SVNRevision(SConfig.Node):
    def __init__(self, scons_env, scons_opts, required=False):
        SConfig.Node.__init__(self, scons_env, scons_opts, required)
        self.checks = [self.extract_revision]
        self.define_name = 'VERSION'
        self.checkout_path = os.getcwd()

        # Will be set after configuration.
        self.revision = 0

    def extract_revision(self):
        svn_path = os.path.join(self.checkout_path, '.svn', 'entries')
        if not os.path.exists(svn_path):
            return [0, '', 'Could not find .svn directory']
        f = file(svn_path, 'r')
	all_lines = f.readlines()
        f.close()

	for l in all_lines:
            ind = l.rfind('revision=')
            if ind != -1:
                self.revision = int(l[ind + 10:l.rfind('"')])
                return True

        self.revision = int(all_lines[3])
        return True
