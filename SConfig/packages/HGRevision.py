import os
import SConfig

class HGRevision(SConfig.Node):
    def __init__(self, scons_env, scons_opts, required=False):
        SConfig.Node.__init__(self, scons_env, scons_opts, required)
        self.checks = [self.extract_revision]
        self.define_name = 'VERSION'

        # Will be set after configuration.
        self.revision = ''

    def extract_revision(self):
        import commands
        result = commands.getstatusoutput('hg identify')
        if result[0]:
            return False
        self.revision = result[1].split()[0].strip()
        return True

    def enable(self, scons_env, old_state=None):
        SConfig.Node.enable(self, scons_env, old_state)
        self.backup_variable(scons_env, 'CPPDEFINES', old_state)
        ver = scons_env['ESCAPE']('"' + self.revision + '"')
        scons_env.AppendUnique(CPPDEFINES=[(self.define_name, ver)])
